/*
 * The cado-nfs dashboard.
 *
 * A hash-routed page with four views, each of which is a pure function
 * of what the api last returned. No framework, no build step, and no
 * third-party code: everything here is part of cado-nfs.
 *
 * Updates are polled rather than streamed. The werkzeug server runs a
 * handful of threads (min(2*ncpu, 4) + 1), and a long-lived
 * server-sent-events connection per open tab would eat them. Polling
 * with If-None-Match costs a 304 when nothing has changed, which is
 * most of the time.
 */

import * as api from './api.js';
import {
    h, append, clear, card, banner, empty, table, pill, figure, tokenGate,
    duration, count, percent, ago,
} from './components.js';
import {dial, stackedBar, legend, barRows, miniBar,
        sparkline} from './charts.js';

/* Poll intervals, in milliseconds. The overview is what people leave
 * open, so it refreshes briskly; the rest is on demand. */
const FAST = 2000;
const SLOW = 10000;

/* ... but not at the computation's expense. cado-nfs.py serves clients
 * from a single thread unless server.threaded is set, so every request
 * this page makes is a request some client is not being served. When
 * the server takes a long time to answer -- which is exactly when it is
 * busy -- wait proportionally longer before asking again. Observed on a
 * 1400-client run: answers took seconds, and polling every two seconds
 * regardless would have been taking a real bite out of the pool. */
const BACKOFF_FACTOR = 4;
const MAX_INTERVAL = 60000;

const WU_SEGMENTS = [
    {key: 'VERIFIED_OK', label: 'verified', colour: '--ok'},
    {key: 'RECEIVED_OK', label: 'received', colour: '--idle'},
    {key: 'ASSIGNED', label: 'assigned', colour: '--accent'},
    {key: 'AVAILABLE', label: 'available', colour: '--faint'},
    {key: 'NEED_RESUBMIT', label: 'to resubmit', colour: '--warn'},
    {key: 'VERIFIED_ERROR', label: 'failed', colour: '--bad'},
    {key: 'RECEIVED_ERROR', label: 'errored', colour: '--bad'},
    {key: 'CANCELLED', label: 'cancelled', colour: '--border'},
];

/* How the pool can be rolled up, outermost first. */
const GROUPINGS = ['domain', 'cluster', 'host'];

const VIEWS = [
    {id: 'overview', label: 'Overview'},
    {id: 'clients', label: 'Clients'},
    {id: 'workunits', label: 'Workunits'},
    {id: 'log', label: 'Log'},
];

const CLIENT_PAGE = 100;

const state = {
    view: 'overview',
    id: null,              /* what the view is about, if anything */
    params: new URLSearchParams(),
    info: null,
    progress: null,
    summary: null,
    clients: null,
    clientsSummary: null,
    workunits: null,
    log: null,
    error: null,
    notice: null,
    lastOk: null,
    busy: 0,
    clientSort: {sort: 'completed', desc: true},
    parameters: null,
    detail: null,          /* the drill-down currently on screen */
};

let timer = null;

/* ---------------- routing ----------------
 *
 * The hash is the whole of the page's state: #view/id?filters. Every
 * filter control and every drill-down therefore produces a link that
 * can be bookmarked, sent to a colleague, or opened in a second tab --
 * and the back button does what a back button should. */

function parseHash() {
    const raw = (location.hash || '#overview').replace(/^#/, '');
    const cut = raw.indexOf('?');
    const path = cut < 0 ? raw : raw.slice(0, cut);
    const segments = path.split('/').filter((s) => s !== '')
        .map(decodeURIComponent);
    return {
        view: segments[0] || 'overview',
        id: segments.length > 1 ? segments.slice(1).join('/') : null,
        params: new URLSearchParams(cut < 0 ? '' : raw.slice(cut + 1)),
    };
}

/* Build a hash. Values that are empty are left out, so that a link
 * carries only the filters that are actually on. */
function link(view, id, params) {
    let hash = '#' + view;
    if (id !== null && id !== undefined && id !== '') {
        hash += '/' + encodeURIComponent(id);
    }
    const query = new URLSearchParams();
    for (const [key, value] of Object.entries(params || {})) {
        if (value !== null && value !== undefined && value !== ''
            && value !== false) {
            query.set(key, value);
        }
    }
    const q = query.toString();
    return hash + (q ? '?' + q : '');
}

/* A link that keeps the filters currently in force and changes some. */
function here(changes) {
    const params = {};
    for (const [key, value] of state.params.entries()) params[key] = value;
    Object.assign(params, changes);
    /* Any change to a filter starts the paging over: page 4 of the old
     * answer has nothing to do with the new one. */
    if (!('offset' in changes)) delete params.offset;
    return link(state.view, state.id, params);
}

function go(hash) {
    if (location.hash === hash) {
        refresh(true);
    } else {
        location.hash = hash;
    }
}

function param(name, fallback = '') {
    const value = state.params.get(name);
    return value === null || value === undefined ? fallback : value;
}

function intParam(name, fallback) {
    const value = Number(state.params.get(name));
    return isNaN(value) || !state.params.get(name) ? fallback : value;
}

/* ---------------- chrome ---------------- */

/* Which roll-up to show when nobody has said. The useful default is
 * the outermost one that both discriminates and compresses: grouping
 * by domain says nothing when every client is in the same domain, and
 * grouping by machine says nothing when every machine runs exactly one
 * client -- that is the flat list with columns taken away. */
function defaultGrouping() {
    const summary = state.clientsSummary || {};
    const census = summary.groupings;
    if (!census) return 'cluster';
    for (const kind of GROUPINGS) {
        const groups = census[kind] || 0;
        if (groups > 1 && groups < summary.total) return kind;
    }
    return '';
}

function showBusy() {
    const bar = document.getElementById('busy');
    if (bar) bar.className = 'busy' + (state.busy > 0 ? ' on' : '');
}

function heartbeat() {
    const dot = document.getElementById('heartbeat');
    if (!dot) return;
    const stale = state.lastOk !== null
        && (Date.now() - state.lastOk) > 4 * SLOW;
    dot.className = 'heartbeat' + (stale ? ' stale' : '');
    if (!stale) {
        dot.classList.add('beat');
        setTimeout(() => dot.classList.remove('beat'), 200);
    }
    dot.title = (state.lastOk
        ? 'last successful update ' + duration((Date.now() - state.lastOk)
                                               / 1000) + ' ago'
        : 'no update yet')
        + '; last round trip ' + Math.round(api.lastRoundTripMs)
        + ' ms, polling every ' + Math.round(pollInterval() / 1000) + ' s';
}

function renderTopbar() {
    const info = state.info || {};
    const progress = state.progress || {};
    const bar = document.getElementById('topbar');
    clear(bar);

    const meta = [];
    if (info.starttime) {
        meta.push(h('span', {}, 'running for ',
                    h('strong', {},
                      duration(api.now() - Number(info.starttime)))));
    }
    const kind = info.computation_desc || info.computation;
    if (kind) meta.push(h('span', {}, kind));
    if (progress.finished) {
        meta.push(pill('finished', 'working'));
    } else if (info.serving_workunits === false) {
        meta.push(pill('not serving workunits', 'stale'));
    }

    append(bar, [
        h('div', {class: 'brand'}, 'cado-nfs ',
          h('span', {class: 'name'}, info.name || '')),
        h('div', {class: 'meta'}, meta),
        h('div', {class: 'spacer'}),
        h('span', {id: 'busy', class: 'busy'}),
        h('span', {id: 'heartbeat', class: 'heartbeat'}),
        h('button', {
            class: 'small',
            title: 'forget the api token in this tab',
            onclick: () => { api.forgetToken(); location.reload(); },
        }, 'Lock'),
    ]);
    heartbeat();
}

/* A drill-down still belongs to the tab it was reached from. */
const TAB_OF = {group: 'clients', stage: 'overview'};

function renderTabs() {
    const nav = document.getElementById('tabs');
    clear(nav);
    const active = TAB_OF[state.view] || state.view;
    for (const view of VIEWS) {
        nav.appendChild(h('a', {
            href: '#' + view.id,
            class: view.id === active ? 'active' : '',
        }, view.label));
    }
}

/* ---------------- overview ---------------- */

function currentTask() {
    const tasks = (state.progress && state.progress.tasks) || [];
    return tasks.find((t) => t.phase === 'running') || null;
}

function highlights(task) {
    const bits = [];
    const hl = (task && task.highlights) || {};
    if (hl.rels_found !== undefined && hl.rels_wanted !== undefined) {
        bits.push(figure(count(hl.rels_found) + ' / ' + count(hl.rels_wanted),
                         'relations'));
    }
    if (hl.qnext !== undefined) {
        bits.push(figure(count(hl.qnext), 'special-q reached'));
    }
    if (hl.adnext !== undefined) {
        bits.push(figure(count(hl.adnext), 'leading coefficient'));
    }
    if (task && task.wu_submitted !== undefined) {
        bits.push(figure(
            h('a', {href: link('workunits', null, {task: task.name})},
              count(task.wu_received || 0) + ' / '
              + count(task.wu_submitted)),
            'workunits back'));
    }
    if (task && task.wu_failed) {
        bits.push(figure(
            h('a', {href: link('workunits', null,
                               {task: task.name,
                                status: 'VERIFIED_ERROR'}),
                    style: 'color:var(--bad)'}, count(task.wu_failed)),
            'failed'));
    }
    return bits;
}

function currentCard() {
    const task = currentTask();
    if (!task) {
        const finished = state.progress && state.progress.finished;
        return card('Current phase',
                    empty(finished
                          ? 'The computation has finished.'
                          : 'No task is running yet.'));
    }
    const achievement = task.achievement;
    return card('Current phase',
                h('div', {class: 'current'},
                  h('div', {class: 'dial'},
                    dial(achievement,
                         achievement === undefined ? '–'
                             : percent(achievement, 1),
                         'complete')),
                  h('div', {class: 'detail'},
                    h('div', {class: 'title'},
                      h('a', {href: link('stage', task.name)},
                        task.title || task.name)),
                    h('div', {class: 'sub'},
                      task.eta ? 'estimated finish ' + task.eta
                          : 'no estimate yet',
                      state.progress.current_started
                          ? ' · started ' + ago(
                              state.progress.current_started, api.now())
                              + ' ago'
                          : ''),
                    h('div', {class: 'figures'}, highlights(task)))));
}

/* What the right-hand end of a pipeline row says, and why. */
const PHASE_NOTE = {
    disabled: ['not run', 'tasks.<name>.run is false for this run, and'
               + ' a task turned off that way stops the run rather'
               + ' than being stepped over'],
    unreachable: ['not reached', 'the run stops at a disabled task'
                  + ' before this one'],
    pending: ['', 'not started yet'],
};

function pipelineCard() {
    const tasks = (state.progress && state.progress.tasks) || [];
    if (!tasks.length) return card('Pipeline', empty('Nothing published yet.'));
    const list = h('ol', {class: 'pipeline'});
    for (const task of tasks) {
        let right = '';
        if (task.phase === 'done') {
            right = 'done';
        } else if (task.phase === 'running'
                   && task.achievement !== undefined) {
            right = percent(task.achievement, 1);
        } else if (PHASE_NOTE[task.phase]) {
            right = PHASE_NOTE[task.phase][0];
        }
        list.appendChild(h('li', {
            class: task.phase,
            title: (PHASE_NOTE[task.phase] || [])[1] || '',
        },
            h('span', {class: 'dot'}),
            h('a', {class: 'label', href: link('stage', task.name)},
              task.title || task.name),
            h('span', {class: 'right'}, right)));
    }
    return card('Pipeline', list);
}

function workunitSegments() {
    const counts = (state.summary && state.summary.counts) || {};
    return WU_SEGMENTS.map((s) => ({
        label: s.label,
        colour: s.colour,
        value: counts[s.key] || 0,
        hideEmpty: true,
        /* Every tally on this page is a question about the rows behind
         * it, so every tally is a link to them. */
        href: link('workunits', null, {status: s.key}),
    }));
}

function workunitsCard() {
    const segments = workunitSegments();
    const total = (state.summary && state.summary.total) || 0;
    return card('Workunits',
                h('div', {class: 'figures'},
                  figure(h('a', {href: link('workunits')}, count(total)),
                         'total'),
                  figure(h('a', {href: link('workunits', null,
                                            {status: 'ASSIGNED'})},
                           count((state.summary || {}).outstanding || 0)),
                         'outstanding')),
                h('div', {style: 'margin-top:14px'},
                  stackedBar(segments)),
                legend(segments));
}

function clientsSummaryCard() {
    /* Fed by /api/v1/clients?summary=1, never by the full list: on a
     * pool of any size that list is far too big to poll. */
    const payload = state.clientsSummary || {};
    const rows = payload.top || [];
    if (!payload.total) {
        return card('Clients', empty('No client has asked for work yet.'));
    }
    const counts = payload.counts || {};
    const order = ['working', 'idle', 'stale', 'gone', 'unknown'];
    const chips = order.filter((k) => counts[k])
        .map((k) => h('a', {href: link('clients', null,
                                       {state: k, group_by: 'flat'}),
                            class: 'pilllink'},
                      pill(counts[k] + ' ' + k, k)));

    const top = rows.map((c) => ({
        label: c.clientid,
        value: c.completed,
        href: link('clients', c.clientid),
        colour: c.state === 'gone' ? '--bad'
            : c.state === 'stale' ? '--warn' : '--accent',
    }));

    return card('Clients',
                h('div', {}, chips, ' ',
                  h('a', {href: link('clients'), class: 'faint'},
                    payload.total + ' in all')),
                h('div', {style: 'margin-top:12px'},
                  barRows(top, {format: count})),
                payload.truncated
                    ? h('div', {class: 'faint',
                                style: 'margin-top:8px;font-size:12px'},
                        'and ' + payload.truncated + ' more \u2014 ',
                        h('a', {href: link('clients')}, 'see all'))
                    : null);
}

function strandedCard() {
    const summary = state.clientsSummary || {};
    const counts = summary.counts || {};
    const quiet = (counts.stale || 0) + (counts.gone || 0);
    if (!quiet) return null;
    const rows = ((state.clients || {}).clients || [])
        .filter((c) => c.in_flight
                && (c.state === 'stale' || c.state === 'gone'));
    const quietList = link('clients', null,
                           {state: 'stale,gone', group_by: 'flat'});
    if (!rows.length) {
        /* Naming them needs the full list; point at exactly those
         * clients rather than at the whole pool, which on a large run
         * is thousands of rows one would then have to sift. */
        return h('section', {class: 'card span2'},
                 banner('warn',
                        h('strong', {}, quiet + ' client'
                          + (quiet === 1 ? '' : 's')),
                        ' have gone quiet. ',
                        h('a', {href: quietList},
                          'List those ' + quiet),
                        ' to see what they are holding.'));
    }
    const held = rows.reduce((s, c) => s + c.in_flight, 0);
    return h('section', {class: 'card span2'},
             banner('warn',
                    h('strong', {},
                      held + ' workunit' + (held === 1 ? '' : 's')),
                    ' held by ' + rows.length + ' client'
                    + (rows.length === 1 ? '' : 's')
                    + ' we have not heard from recently.'),
             table([
                 {key: 'clientid', label: 'client', mono: true,
                  render: (r) => h('a', {href: link('clients', r.clientid)},
                                   r.clientid)},
                 {key: 'state', label: 'state',
                  render: (r) => statePill(r)},
                 {key: 'in_flight', label: 'holding', num: true},
                 {key: 'typical_turnaround', label: 'usual pace',
                  num: true,
                  render: (r) => r.typical_turnaround === null
                      ? h('span', {class: 'faint'}, '–')
                      : duration(r.typical_turnaround)},
                 {key: 'last_seen', label: 'last seen', num: true,
                  render: (r) => ago(r.last_seen, api.now())},
                 {key: 'act', label: '', sort: false,
                  render: (r) => h('button', {
                      class: 'small',
                      onclick: (e) => reclaim(e.target, r.clientid),
                  }, 'Reclaim')},
             ], rows),
             h('div', {class: 'controls', style: 'margin-top:12px'},
               h('button', {
                   onclick: (e) => stateReclaim(e.target, 'stale', 'gone'),
               }, 'Reclaim everything these ' + rows.length
                  + ' are holding'),
               h('a', {href: quietList, class: 'button small'},
                 'List all ' + quiet + ' quiet clients')),
             h('div', {class: 'controls', style: 'margin-top:8px'},
               h('button', {
                   onclick: (e) => bulkReclaim(e.target),
               }, 'Reclaim everything assigned longer ago than'),
               h('input', {type: 'number', id: 'bulk-minutes', value: 60,
                           min: 1, style: 'width:80px'}),
               h('span', {class: 'muted'}, 'minutes')));
}

/* The two ceilings whose only effect is to abort the computation. They
 * live on the overview because the moment you want them is the moment
 * a long run is about to be lost. */
function ceilingsCard() {
    const tunables = (state.parameters || {}).parameters;
    if (!tunables) return null;
    const rows = Object.keys(tunables).sort().map((name) => {
        const p = tunables[name];
        const near = p.counter_value !== null
            && p.counter_value !== undefined
            && p.value > 0 && p.counter_value / p.value >= 0.8;
        return h('tr', {},
            h('td', {class: 'mono'}, name),
            h('td', {class: 'num'},
              h('span', {style: near ? 'color:var(--bad);font-weight:600'
                         : ''},
                (p.counter_value === null || p.counter_value === undefined
                 ? '\u2013' : p.counter_value) + ' / ' + p.value),
              p.overridden
                  ? h('span', {class: 'faint'}, ' (raised)') : null,
              /* The counters belong to client-server tasks. While one
               * that keeps none is running, the number shown is the
               * one the last task finished on, and saying so is the
               * difference between a stale figure and a wrong one. */
              p.counter_from && !p.counter_is_current
                  ? h('a', {href: link('stage', p.counter_from),
                            class: 'faint'},
                      ' (' + p.counter_from + ')')
                  : null),
            h('td', {class: 'num'},
              h('button', {
                  class: 'small',
                  title: 'Double this ceiling. It only decides when the'
                      + ' computation gives up, so raising it cannot'
                      + ' change what is computed \u2014 and the change'
                      + ' is written to a new parameters snapshot.',
                  onclick: (e) => raiseCeiling(e.target, name,
                                               2 * p.value),
              }, 'double')));
    });
    return card('Give-up thresholds',
        h('div', {class: 'tablewrap'}, h('table', {}, h('tbody', {}, rows))),
        h('p', {class: 'faint', style: 'font-size:12px;margin-bottom:0'},
          'These abort the computation when the counter reaches them. ',
          'Raising one changes nothing about what is computed, and is ',
          'recorded in a fresh parameters snapshot so the run stays ',
          'reproducible.'));
}

function raiseCeiling(button, name, value) {
    return guard(button, async () => {
        const r = await api.setParameter(name, value);
        return {marked: [name], skipped: [], message: r.message};
    });
}

function overview() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const stranded = strandedCard();
    append(main, [
        h('div', {class: 'grid wide'}, currentCard(), workunitsCard()),
        h('div', {style: 'height:16px'}),
        h('div', {class: 'grid wide'}, pipelineCard(), clientsSummaryCard()),
        stranded ? h('div', {style: 'height:16px'}) : null,
        stranded ? h('div', {class: 'grid'}, stranded) : null,
        h('div', {style: 'height:16px'}),
        h('div', {class: 'grid wide'}, ceilingsCard() || h('div', {})),
    ]);
}

/* ---------------- clients ---------------- */

/* Roll-up or flat list. On a pool of any size the roll-up is the only
 * readable option, and it is also the only bounded one. */
function groupControl() {
    const current = groupBy() === null ? defaultGrouping() : groupBy();
    return h('label', {}, 'group by',
        h('select', {
            onchange: (e) => go(here({group_by: e.target.value || 'flat'})),
        }, [['', 'nothing (list every client)'],
            ['host', 'machine'],
            ['cluster', 'cluster'],
            ['domain', 'domain']].map(([v, label]) => h('option', {
            value: v, selected: v === current,
        }, label))));
}

/* The liveness filter, which is also where "816 clients have gone
 * quiet" lands. */
function stateControl() {
    const counts = (state.clientsSummary || {}).counts || {};
    const current = param('state');
    const choices = [['', 'any state']].concat(
        ['working', 'idle', 'stale', 'gone', 'unknown']
            .filter((k) => counts[k])
            .map((k) => [k, k + ' (' + counts[k] + ')']));
    if (current && !choices.some(([v]) => v === current)) {
        choices.push([current, current]);
    }
    return h('label', {}, 'state',
        h('select', {
            onchange: (e) => go(here({state: e.target.value})),
        }, choices.map(([v, label]) => h('option', {
            value: v, selected: v === current,
        }, label))));
}

/* Every pill in a roll-up is a way in: the point of a tally is the
 * rows behind it. */
function groupStatePills(kind, group) {
    const states = group.states || {};
    const pills = Object.keys(states).sort().map((k) => h('a', {
        href: link('group', group.key, {group_by: kind, state: k}),
        class: 'pilllink',
    }, pill(states[k] + ' ' + k, k)));
    if (group.completed) {
        pills.push(h('a', {
            href: link('workunits', null,
                       {client: group.key, status: 'VERIFIED_OK'}),
            class: 'pilllink',
            title: 'workunits this ' + kind + ' has had verified',
        }, pill(count(group.completed) + ' ok', 'ok')));
    }
    return pills;
}

function clientGroupsView() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const payload = state.clients || {};
    const groups = payload.groups || [];
    const kind = payload.group_by || defaultGrouping();
    main.appendChild(card('Clients by ' + kind,
        h('div', {class: 'controls'}, groupControl(), stateControl(),
          h('span', {class: 'muted'},
            (payload.total || 0) + ' clients in '
            + (payload.groups_total || 0) + ' ' + kind + 's')),
        groups.length
            ? table([
                {key: 'key', label: kind, mono: true,
                 render: (g) => h('a', {
                     href: link('group', g.key, {group_by: kind}),
                 }, g.key)},
                {key: 'clients', label: 'clients', num: true},
                {key: 'cores', label: 'cores', num: true,
                 render: (g) => g.cores
                     || h('span', {class: 'faint'}, '–')},
                {key: 'in_flight', label: 'in flight', num: true},
                {key: 'completed', label: 'completed', num: true},
                {key: 'failed', label: 'failed', num: true,
                 render: (g) => g.failed
                     ? h('span', {style: 'color:var(--bad)'}, g.failed)
                     : h('span', {class: 'faint'}, '0')},
                {key: 'share', label: 'share', num: true,
                 render: (g) => percent(g.share, 0)},
                {key: 'states', label: 'states', sort: false,
                 render: (g) => groupStatePills(kind, g)},
            ], groups, {state: state.clientSort,
                        onsort: () => clientGroupsView()})
            : empty('Nothing to group yet.')));
}

/* ---------------- one machine, one cluster, one domain ------------ */

function groupView() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const kind = param('group_by', 'cluster');
    const payload = state.clients || {};
    const rows = payload.clients || [];
    const filtered = param('state');

    const held = rows.reduce((s, c) => s + c.in_flight, 0);
    const tallies = {};
    for (const c of rows) tallies[c.state] = (tallies[c.state] || 0) + 1;

    main.appendChild(card(null,
        h('div', {}, backLink(link('clients', null, {group_by: kind}),
                              'all ' + kind + 's')),
        h('div', {class: 'detail', style: 'margin-top:10px'},
          h('div', {class: 'title mono'}, state.id,
            h('span', {class: 'faint'}, ' \u00b7 ' + kind)),
          h('div', {class: 'sub'},
            Object.keys(tallies).sort().map(
                (k) => h('a', {href: here({state: k}), class: 'pilllink'},
                         pill(tallies[k] + ' ' + k, k))),
            filtered
                ? h('a', {href: here({state: ''}), class: 'pilllink'},
                    pill('showing ' + filtered + ' only \u2014 clear',
                         'unknown'))
                : null),
          h('div', {class: 'figures'},
            figure(count(rows.length), 'clients'),
            figure(count(rows.reduce((s, c) => s + c.completed, 0)),
                   'completed'),
            figure(count(rows.reduce((s, c) => s + c.failed, 0)),
                   'failed'),
            figure(count(held), 'in flight'),
            figure(percent(rows.reduce((s, c) => s + c.share, 0), 0),
                   'share of all work'))),
        held
            ? h('div', {class: 'controls', style: 'margin-top:14px'},
                h('button', {
                    title: 'Put back in the pool everything the clients'
                        + ' listed here are holding. Any of them that'
                        + ' is alive has its upload refused as a'
                        + ' duplicate and moves on; nothing is lost.',
                    onclick: (e) => groupReclaim(e.target, kind),
                }, 'Reclaim the ' + held + ' workunit'
                   + (held === 1 ? '' : 's') + ' held here'))
            : null));

    main.appendChild(h('div', {style: 'height:16px'}));
    main.appendChild(card('Clients',
        rows.length ? clientTable(rows) : empty('Nothing here.')));
}

/* The list of clients, used both by the flat view and by a group's
 * own page. */
function clientTable(rows, redraw) {
    const totalDone = rows.reduce((s, c) => s + c.completed, 0) || 1;
    return table([
        {key: 'clientid', label: 'client', mono: true,
         render: (r) => h('a', {href: link('clients', r.clientid)},
                          r.clientid)},
        {key: 'state', label: 'state',
         render: (r) => statePill(r)},
        {key: 'typical_turnaround', label: 'usual pace', num: true,
         render: (r) => r.typical_turnaround === null
             ? h('span', {class: 'faint',
                          title: 'not enough recent workunits'
                                 + ' from this client yet'}, '\u2013')
             : h('span', {title: 'median of '
                                 + r.turnaround_samples
                                 + ' recent workunits'},
                 duration(r.typical_turnaround))},
        {key: 'in_flight', label: 'in flight', num: true,
         render: (r) => r.in_flight
             ? h('a', {href: link('workunits', null,
                                  {client: r.clientid,
                                   status: 'ASSIGNED'})}, r.in_flight)
             : h('span', {class: 'faint'}, '0')},
        {key: 'completed', label: 'completed', num: true},
        {key: 'failed', label: 'failed', num: true,
         render: (r) => r.failed
             ? h('a', {href: link('workunits', null,
                                  {client: r.clientid,
                                   status: 'VERIFIED_ERROR'}),
                       style: 'color:var(--bad)'}, r.failed)
             : h('span', {class: 'faint'}, '0')},
        {key: 'share', label: 'share', num: true,
         render: (r) => h('span', {},
                          miniBar(r.completed / totalDone),
                          ' ' + percent(r.share, 0))},
        {key: 'last_seen', label: 'last seen', num: true,
         render: (r) => ago(r.last_seen, api.now())},
        {key: 'oldest_assignment', label: 'oldest task', num: true,
         render: (r) => ago(r.oldest_assignment, api.now())},
        {key: 'act', label: '', sort: false,
         render: (r) => r.in_flight
             ? h('button', {
                 class: 'small',
                 title: 'Put this client\'s workunits back in the'
                     + ' pool. If it turns out to be alive, its'
                     + ' upload is refused as a duplicate and it'
                     + ' moves on; nothing is lost.',
                 onclick: (e) => reclaim(e.target, r.clientid),
             }, 'Reclaim')
             : null},
    ], rows, {state: state.clientSort, onsort: redraw || (() => render())});
}

function pager(total, shown) {
    const offset = intParam('offset', 0);
    if (total <= shown && offset === 0) return null;
    return h('div', {class: 'controls', style: 'margin-bottom:0'},
        h('a', {
            class: 'button small' + (offset <= 0 ? ' disabled' : ''),
            href: offset <= 0 ? null
                : here({offset: Math.max(0, offset - CLIENT_PAGE)}),
        }, '\u2190 previous'),
        h('span', {class: 'muted'},
          (offset + 1) + '\u2013' + (offset + shown) + ' of ' + total),
        h('a', {
            class: 'button small'
                + (offset + shown >= total ? ' disabled' : ''),
            href: offset + shown >= total ? null
                : here({offset: offset + CLIENT_PAGE}),
        }, 'next \u2192'));
}

function clientsView() {
    if (groupBy() === null ? defaultGrouping() : groupBy()) {
        return clientGroupsView();
    }
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const payload = state.clients || {};
    const rows = payload.clients || [];
    const filtered = param('state');
    if (!rows.length) {
        main.appendChild(card('Clients',
            h('div', {class: 'controls'}, groupControl(), stateControl()),
            empty(filtered
                  ? 'No client is ' + filtered + '.'
                  : 'No client has asked for work yet.')));
        return;
    }
    const total = payload.total || rows.length;
    const page = pager(total, rows.length);
    main.appendChild(card(
        'Clients',
        h('div', {class: 'controls'}, groupControl(), stateControl(),
          filtered
              ? h('span', {class: 'muted'},
                  total + ' of ' + (payload.pool_total || total)
                  + ' clients')
              : null,
          filtered && rows.some((c) => c.in_flight)
              ? h('button', {
                  class: 'small',
                  onclick: (e) => stateReclaim(e.target, filtered),
              }, 'Reclaim what they are holding')
              : null),
        h('p', {class: 'muted', style: 'margin-top:-6px'},
          'How long a client may be silent before it counts as stale is ',
          'judged per client, from how long its own workunits have ',
          'recently been taking \u2014 hover a state to see the ',
          'reasoning. It is never longer than tasks.wutimeout (',
          duration(payload.wutimeout), '), which is when its work gets ',
          'reassigned anyway, and while few workunits back the ',
          'estimate it is held down to a couple of ',
          'tasks.wutimeoutcheck intervals.'),
        page,
        clientTable(rows, () => clientsView()),
        page ? h('div', {class: 'faint',
                         style: 'margin-top:10px;font-size:12px'},
                 'Sorting applies to this page only: the server orders'
                 + ' clients by what they have contributed, and the'
                 + ' dashboard asks for one page at a time so that a'
                 + ' large pool does not cost the computation'
                 + ' bandwidth it needs for workunits.')
             : null));
}

/* ---------------- drill-downs ---------------- */

function backLink(href, label) {
    return h('a', {href, style: 'font-size:13px'}, '\u2190 ' + label);
}

function timeline(entries) {
    const rows = entries.filter((e) => e[1])
        .map(([label, stamp]) => h('tr', {},
            h('td', {class: 'faint'}, label),
            h('td', {class: 'num'}, ago(stamp, api.now()) + ' ago')));
    if (!rows.length) return null;
    return h('div', {class: 'tablewrap'},
             h('table', {}, h('tbody', {}, rows)));
}

function clientDetailView() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const c = state.detail;
    if (!c) {
        main.appendChild(card('Client', empty('Loading\u2026')));
        return;
    }

    const pace = c.typical_turnaround === null
        ? h('span', {class: 'faint'}, 'not enough workunits yet')
        : h('span', {}, duration(c.typical_turnaround),
            h('span', {class: 'faint'},
              ' (median of ' + c.turnaround_samples + ')'));

    main.appendChild(card(null,
        h('div', {}, backLink(link('clients'), 'all clients')),
        h('div', {class: 'current', style: 'margin-top:10px'},
          h('div', {class: 'detail'},
            h('div', {class: 'title mono'}, c.clientid, ' ', statePill(c)),
            h('div', {class: 'sub'}, livenessTitle(c)),
            h('div', {class: 'sub faint'},
              /* Where it sits is a place one wants to go: several
               * clients share a machine, and machines come and go a
               * rack at a time. */
              c.host
                  ? h('a', {href: link('group', c.host,
                                       {group_by: 'host'})},
                      c.fqdn || c.host)
                  : (c.fqdn || '\u2013'),
              c.cluster
                  ? h('span', {}, ' \u00b7 cluster ',
                      h('a', {href: link('group', c.cluster,
                                         {group_by: 'cluster'})},
                        c.cluster))
                  : '',
              c.domain
                  ? h('span', {}, ' \u00b7 ',
                      h('a', {href: link('group', c.domain,
                                         {group_by: 'domain'})},
                        c.domain))
                  : '',
              c.cores ? ' \u00b7 ' + c.cores + ' cores' : '',
              c.platform ? ' \u00b7 ' + c.platform : '',
              Object.keys(c.overrides || {}).length
                  ? ' \u00b7 overrides ' + Object.entries(c.overrides)
                      .map(([k, v]) => k + '=' + v).join(' ')
                  : '',
              c.self_reported ? '' : ' (inferred from its name)'),
            h('div', {class: 'figures'},
              figure(count(c.completed), 'completed'),
              figure(count(c.failed), 'failed'),
              figure(count(c.in_flight), 'in flight'),
              figure(percent(c.share, 0), 'share of all work'),
              figure(pace, 'usual pace'),
              figure(ago(c.last_seen, api.now()), 'last seen'))),
          c.turnarounds && c.turnarounds.length > 1
              ? h('div', {},
                  h('div', {class: 'label',
                            style: 'text-align:right'},
                    'recent turnaround'),
                  sparkline(c.turnarounds.slice().reverse()))
              : null),
        c.in_flight
            ? h('div', {style: 'margin-top:14px'},
                h('button', {
                    onclick: (e) => reclaim(e.target, c.clientid),
                }, 'Reclaim its ' + c.in_flight + ' workunit'
                   + (c.in_flight === 1 ? '' : 's')))
            : null));

    const columns = [
        {key: 'wuid', label: 'workunit', mono: true,
         render: (w) => h('a', {href: '#workunits/'
                                + encodeURIComponent(w.wuid)}, w.wuid)},
        {key: 'status_name', label: 'status',
         render: (w) => statusPill(w.status_name)},
        {key: 'attempt', label: 'try', num: true},
        {key: 'duration', label: 'took', num: true,
         render: (w) => duration(w.duration)},
        {key: 'timeresult', label: 'returned', num: true,
         render: (w) => ago(w.timeresult || w.timeassigned, api.now())},
    ];

    main.appendChild(h('div', {style: 'height:16px'}));
    main.appendChild(card('Holding now',
        (c.in_flight_workunits || []).length
            ? table(columns, c.in_flight_workunits, {state: {}})
            : empty('Nothing.')));
    main.appendChild(h('div', {style: 'height:16px'}));
    main.appendChild(card('Recently returned',
        (c.recent_workunits || []).length
            ? table(columns, c.recent_workunits, {state: {}})
            : empty('Nothing yet.')));
}

function workunitDetailView() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const w = state.detail;
    if (!w) {
        main.appendChild(card('Workunit', empty('Loading\u2026')));
        return;
    }

    const who = (id) => id
        ? h('a', {href: '#clients/' + encodeURIComponent(id),
                  class: 'mono'}, id)
        : h('span', {class: 'faint'}, '\u2013');

    main.appendChild(card(null,
        h('div', {}, backLink(link('workunits'), 'all workunits')),
        h('div', {style: 'margin-top:10px'},
          h('div', {class: 'title mono'}, w.wuid, ' ',
            statusPill(w.status_name)),
          h('div', {class: 'sub'},
            'task ', h('strong', {}, w.task || '?'),
            ', range ', h('strong', {}, w.identifier || '?'),
            ', attempt ', h('strong', {}, String(w.attempt)))),
        h('div', {class: 'figures'},
          figure(duration(w.duration), 'time held'),
          figure(who(w.assignedclient), 'assigned to'),
          figure(who(w.resultclient), 'returned by'),
          w.errorcode ? figure(w.errorcode, 'exit code') : null)));

    main.appendChild(h('div', {style: 'height:16px'}));
    main.appendChild(h('div', {class: 'grid wide'},
        card('Attempts', (w.attempts_all || []).length
            ? table([
                {key: 'attempt', label: 'try', num: true},
                {key: 'wuid', label: 'workunit', mono: true,
                 render: (a) => a.wuid === w.wuid
                     ? h('strong', {}, a.wuid)
                     : h('a', {href: '#workunits/'
                                     + encodeURIComponent(a.wuid)},
                         a.wuid)},
                {key: 'status_name', label: 'status',
                 render: (a) => statusPill(a.status_name)},
                {key: 'resultclient', label: 'client',
                 render: (a) => who(a.resultclient || a.assignedclient)},
                {key: 'duration', label: 'took', num: true,
                 render: (a) => duration(a.duration)},
            ], w.attempts_all, {state: {}})
            : empty('Just this one.')),
        card('Timeline', timeline([
            ['created', w.timecreated],
            ['assigned', w.timeassigned],
            ['returned', w.timeresult],
            ['verified', w.timeverified],
        ]) || empty('Not started.'))));

    const commands = ((w.workunit || {}).commands) || [];
    if (commands.length) {
        main.appendChild(h('div', {style: 'height:16px'}));
        main.appendChild(card('Commands',
            h('pre', {class: 'log'}, commands.join('\n'))));
    }

    for (const entry of w.output || []) {
        const pre = h('pre', {class: 'log'},
                      (entry.lines || []).join('\n'));
        main.appendChild(h('div', {style: 'height:16px'}));
        main.appendChild(card(entry.type + ' \u2014 ' + entry.filename,
                              pre));
        /* This is a tail, and whatever went wrong is at the end of it. */
        pre.scrollTop = pre.scrollHeight;
    }

    if ((w.files || []).length) {
        main.appendChild(h('div', {style: 'height:16px'}));
        main.appendChild(card('Files', table([
            {key: 'type', label: 'type'},
            {key: 'filename', label: 'name', mono: true},
            {key: 'path', label: 'path', mono: true},
        ], w.files, {state: {}})));
    }
}

/* ---------------- one stage of the pipeline ---------------- */

function stageView() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const tasks = (state.progress && state.progress.tasks) || [];
    const task = tasks.find((t) => t.name === state.id);
    if (!task) {
        main.appendChild(card('Stage',
            tasks.length
                ? banner('warn', 'There is no stage called ',
                         h('code', {}, state.id), ' in this pipeline.')
                : empty('Loading\u2026')));
        return;
    }

    const note = PHASE_NOTE[task.phase];
    const submitted = task.wu_submitted;
    const back = (task.wu_received || 0);

    main.appendChild(card(null,
        h('div', {}, backLink(link('overview'), 'overview')),
        h('div', {class: 'current', style: 'margin-top:10px'},
          task.achievement === undefined ? null
              : h('div', {class: 'dial'},
                  dial(task.achievement, percent(task.achievement, 1),
                       'complete')),
          h('div', {class: 'detail'},
            h('div', {class: 'title'}, task.title || task.name, ' ',
              pill(task.phase, task.phase === 'running' ? 'working'
                   : task.phase === 'done' ? 'ok'
                   : task.phase === 'disabled' ? 'gone' : 'unknown')),
            h('div', {class: 'sub'},
              task.phase === 'running'
                  ? (task.eta ? 'estimated finish ' + task.eta
                      : 'no estimate yet')
                  : (note ? note[1] : '')),
            h('div', {class: 'figures'}, highlights(task))))));

    if (submitted !== undefined) {
        main.appendChild(h('div', {style: 'height:16px'}));
        main.appendChild(card('Workunits',
            h('div', {class: 'figures'},
              figure(h('a', {href: link('workunits', null,
                                        {task: task.name})},
                       count(submitted)), 'submitted'),
              figure(count(back), 'back'),
              figure(count(task.wu_timedout || 0), 'timed out'),
              figure(task.wu_failed
                     ? h('a', {href: link('workunits', null,
                                          {task: task.name,
                                           status: 'VERIFIED_ERROR'}),
                               style: 'color:var(--bad)'},
                         count(task.wu_failed))
                     : count(0), 'failed'),
              task.wu_range_received === undefined ? null
                  : figure(count(task.wu_range_received), 'range done')),
            submitted
                ? h('div', {style: 'margin-top:14px'},
                    stackedBar([
                        {label: 'back', value: back, colour: '--ok'},
                        {label: 'out', value: Math.max(0, submitted - back),
                         colour: '--accent'},
                    ]))
                : null));
    }

    const times = task.times || {};
    const names = Object.keys(times).sort();
    if (names.length || (task.stats || []).length) {
        main.appendChild(h('div', {style: 'height:16px'}));
        main.appendChild(card('Time and statistics',
            names.length
                ? table([
                    {key: 'what', label: 'program'},
                    {key: 'seconds', label: 'time', num: true,
                     render: (r) => duration(r.seconds)},
                ], names.map((k) => ({
                    what: k.replace(/^(cpu|real)time_/,
                                    (m, p) => p + ' time, '),
                    seconds: times[k],
                })), {state: {}})
                : null,
            (task.stats || []).length
                ? h('pre', {class: 'log'}, task.stats.join('\n'))
                : null));
    }
}

/* ---------------- workunits ---------------- */

const STATUSES = ['', 'AVAILABLE', 'ASSIGNED', 'NEED_RESUBMIT',
                  'RECEIVED_OK', 'RECEIVED_ERROR', 'VERIFIED_OK',
                  'VERIFIED_ERROR', 'CANCELLED'];

function workunitsView() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const payload = state.workunits || {};
    const rows = payload.workunits || [];
    const tasks = ((state.progress || {}).tasks || []).map((t) => t.name);
    const offset = intParam('offset', 0);
    const limit = intParam('limit', 50);

    const filters = h('div', {class: 'controls'},
        h('label', {}, 'status',
          h('select', {
              onchange: (e) => go(here({status: e.target.value})),
          }, STATUSES.map((s) => h('option', {
              value: s, selected: s === param('status'),
          }, s || 'any')))),
        h('label', {}, 'task',
          h('select', {
              onchange: (e) => go(here({task: e.target.value})),
          }, [''].concat(tasks).map((t) => h('option', {
              value: t, selected: t === param('task'),
          }, t || 'any')))),
        h('label', {}, 'client',
          h('input', {
              type: 'text', value: param('client'),
              placeholder: 'any part of its name',
              title: 'Matches anywhere in the client id, so a machine'
                  + ' name finds every process on it.',
              onchange: (e) => go(here({client: e.target.value.trim()})),
          })),
        h('label', {}, 'rows',
          h('select', {
              onchange: (e) => go(here({limit: e.target.value})),
          }, [25, 50, 100, 200].map((n) => h('option', {
              value: n, selected: n === limit,
          }, n)))),
        /* The page is one page of an unknown number: the api answers
         * how many rows it returned, not how many exist. */
        h('a', {
            class: 'button small' + (offset <= 0 ? ' disabled' : ''),
            href: offset <= 0 ? null
                : here({offset: Math.max(0, offset - limit)}),
        }, '\u2190'),
        h('span', {class: 'muted'},
          rows.length
              ? (offset + 1) + '\u2013' + (offset + rows.length)
              : 'nothing here'),
        h('a', {
            class: 'button small' + (rows.length < limit ? ' disabled' : ''),
            href: rows.length < limit ? null
                : here({offset: offset + limit}),
        }, '\u2192'));

    main.appendChild(card('Workunits', filters,
        rows.length
            ? table([
                {key: 'wuid', label: 'workunit', mono: true,
                 render: (r) => h('a', {href: link('workunits', r.wuid)},
                                  r.wuid)},
                {key: 'status_name', label: 'status',
                 render: (r) => h('a', {href: here({status: r.status_name}),
                                        class: 'pilllink'},
                                  statusPill(r.status_name))},
                {key: 'task', label: 'task',
                 render: (r) => r.task
                     ? h('a', {href: here({task: r.task})},
                         r.task)
                     : h('span', {class: 'faint'}, '\u2013')},
                {key: 'attempt', label: 'try', num: true},
                {key: 'assignedclient', label: 'client',
                 render: (r) => {
                     const id = r.assignedclient || r.resultclient;
                     return id
                         ? h('a', {href: link('clients', id)}, id)
                         : h('span', {class: 'faint'}, '\u2013');
                 }},
                {key: 'timeassigned', label: 'assigned', num: true,
                 render: (r) => ago(r.timeassigned, api.now())},
                {key: 'timeresult', label: 'returned', num: true,
                 render: (r) => ago(r.timeresult, api.now())},
                {key: 'act', label: '', sort: false,
                 render: (r) => r.status_name === 'ASSIGNED'
                     ? h('button', {
                         class: 'small',
                         onclick: (e) => resubmit(e.target, r.wuid),
                     }, 'Resubmit')
                     : null},
            ], rows, {state: {}, onsort: () => workunitsView()})
            : empty('No workunit matches these filters.')));
}

/* Why a client is in the state it is in. The threshold is per client
 * and derived from its own history, so asserting "stale" without
 * saying what that was measured against would be unhelpful. */
function livenessTitle(client) {
    const seen = client.last_seen === null || client.last_seen === undefined
        ? 'never heard from'
        : 'last heard from ' + ago(client.last_seen, api.now()) + ' ago';
    if (client.liveness_basis === 'turnaround') {
        return seen + '; usually returns a workunit in '
            + duration(client.typical_turnaround) + ' (median of '
            + client.turnaround_samples + '), so counted stale after '
            + duration(client.stale_after);
    }
    return seen + '; too few recent workunits to know its usual pace,'
        + ' so counted stale after tasks.wutimeout ('
        + duration(client.stale_after) + ')';
}

function statePill(client) {
    const el = pill(client.state, client.state);
    el.title = livenessTitle(client);
    return el;
}

function statusPill(name) {
    const kind = name === 'VERIFIED_OK' || name === 'RECEIVED_OK' ? 'working'
        : name === 'ASSIGNED' ? 'accent'
        : name === 'NEED_RESUBMIT' ? 'stale'
        : name && name.endsWith('_ERROR') ? 'gone'
        : 'unknown';
    return pill(name, kind);
}

/* ---------------- log ---------------- */

function logView() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const payload = state.log || {};
    const pre = h('pre', {class: 'log'});
    for (const line of payload.lines || []) {
        pre.appendChild(h('span', {class: logClass(line)}, line + '\n'));
    }
    main.appendChild(card('Log',
        h('div', {class: 'controls'},
          h('label', {}, 'lines',
            h('select', {
                onchange: (e) => go(here({tail: e.target.value})),
            }, [100, 200, 500, 1000, 2000].map((n) => h('option', {
                value: n, selected: n === intParam('tail', 200),
            }, n)))),
          h('span', {class: 'faint mono'}, payload.path || '')),
        (payload.lines || []).length ? pre : empty('The log is empty.')));
    pre.scrollTop = pre.scrollHeight;
}

function logClass(line) {
    /* cadologger.FileFormatter writes
     *     PID<pid> <asctime> <Level>:<logger name>: message
     * with the level title-cased, so match case-insensitively and only
     * in the prefix, where the level actually lives -- otherwise the
     * word "error" inside a message would colour the whole line. */
    const m = /\b(critical|fatal|error|warning|debug)\s*:/i
        .exec(line.slice(0, 64));
    return m ? 'lvl-' + m[1].toLowerCase() : '';
}

/* ---------------- actions ---------------- */

function notify(kind, message) {
    state.notice = banner(kind, message);
    render();
    setTimeout(() => {
        if (state.notice) { state.notice = null; render(); }
    }, 8000);
}

function describeAction(result) {
    let message = result.message || '';
    if (result.skipped && result.skipped.length) {
        message += ' ' + result.skipped.length + ' skipped: '
            + result.skipped[0].reason;
    }
    return message;
}

async function guard(button, work) {
    if (button) button.disabled = true;
    try {
        const result = await work();
        notify(result.marked && result.marked.length ? 'info' : 'warn',
               describeAction(result));
        await refresh(true);
    } catch (e) {
        notify('error', e.message);
    } finally {
        if (button) button.disabled = false;
    }
}

const reclaim = (button, clientid) =>
    guard(button, () => api.reclaimClient(clientid));

const resubmit = (button, wuid) =>
    guard(button, () => api.resubmitWorkunit(wuid));

/* Everything held by the clients of one machine, cluster or domain. */
const groupReclaim = (button, kind) =>
    guard(button, () => api.reclaimClients({
        group_by: kind, group: state.id,
        states: param('state') ? [param('state')] : undefined,
    }));

/* Everything held by the clients in one liveness state -- which is
 * what "816 clients have gone quiet" is really asking about. */
const stateReclaim = (button, ...states) =>
    guard(button, () => api.reclaimClients({states}));

function bulkReclaim(button) {
    const field = document.getElementById('bulk-minutes');
    const minutes = Math.max(1, Number(field && field.value) || 60);
    return guard(button, () => api.reclaimOlderThan(minutes * 60));
}

/* ---------------- polling ---------------- */

function render() {
    renderTopbar();
    renderTabs();
    if (state.error) {
        const main = document.getElementById('main');
        clear(main);
        main.appendChild(banner('error', state.error));
        return;
    }
    const kind = detailKind();
    if (kind === 'client') return clientDetailView();
    if (kind === 'workunit') return workunitDetailView();
    if (kind === 'group') return groupView();
    if (kind === 'stage') return stageView();
    ({
        overview,
        clients: clientsView,
        workunits: workunitsView,
        log: logView,
    })[state.view]();
}

async function refresh(immediate = false) {
    const started = Date.now();
    /* A request in flight must be visible. On a loaded server the
     * client view takes seconds, and a page that looks identical while
     * it waits is a page that looks broken. */
    state.busy += 1;
    showBusy();
    try {
        /* Only what the view on screen actually needs. Fetching the
         * full client list on every tick regardless -- which is what
         * this used to do -- is half a megabyte every two seconds on a
         * 1400-client run, taken straight out of the server's capacity
         * to hand out workunits. */
        const wanted = [
            api.info().then((r) => { state.info = r; }),
            api.progress().then((r) => { state.progress = r; }),
        ];
        /* Which roll-up to show depends on how many machines,
         * clusters and domains there are, and that census comes with
         * the summary. Once. Afterwards it is polled alongside
         * everything else. */
        if (state.view === 'clients' && state.id === null
            && groupBy() === null && !state.clientsSummary) {
            state.clientsSummary = await api.clientsSummary();
        }

        const kind = detailKind();
        if (state.view === 'overview' || kind === 'stage') {
            wanted.push(api.summary().then((r) => { state.summary = r; }));
            wanted.push(api.clientsSummary()
                .then((r) => { state.clientsSummary = r; }));
            wanted.push(api.parameters()
                .then((r) => { state.parameters = r; }));
        }
        if (kind === 'client') {
            wanted.push(api.client(state.id)
                .then((r) => { state.detail = r; }));
        } else if (kind === 'workunit') {
            wanted.push(api.workunit(state.id)
                .then((r) => { state.detail = r; }));
        } else if (kind === 'group') {
            wanted.push(api.clients(groupQuery())
                .then((r) => { state.clients = r; }));
        } else if (state.view === 'clients') {
            wanted.push(api.clientsSummary()
                .then((r) => { state.clientsSummary = r; }));
            wanted.push(api.clients(clientQuery())
                .then((r) => { state.clients = r; }));
        }
        if (state.view === 'workunits' && kind === null) {
            wanted.push(api.summary().then((r) => { state.summary = r; }));
            wanted.push(api.workunits(workunitQuery())
                .then((r) => { state.workunits = r; }));
        }
        if (state.view === 'log') {
            wanted.push(api.log(intParam('tail', 200))
                .then((r) => { state.log = r; }));
        }
        await Promise.all(wanted);
        state.error = null;
        state.lastOk = Date.now();
    } catch (e) {
        if (e.status === 401) {
            api.forgetToken();
            location.reload();
            return;
        }
        if (e.status === 404 && detailKind()) {
            /* A drill-down onto something that is not there is not an
             * outage; do not leave the previous one on screen. */
            state.detail = null;
        }
        state.error = e.message;
    }
    api.noteRoundTrip(Date.now() - started);
    state.busy -= 1;
    render();
    showBusy();
    if (immediate) schedule();
}

/* The three queries the api takes, each read straight off the hash so
 * that what is on screen and what is in the address bar cannot
 * disagree. */

function groupBy() {
    const wanted = param('group_by');
    if (GROUPINGS.includes(wanted)) return wanted;
    return wanted === 'flat' ? '' : null;   /* null: not chosen yet */
}

function clientQuery() {
    const grouping = groupBy() === null ? defaultGrouping() : groupBy();
    if (grouping) return {group_by: grouping, state: param('state')};
    return {limit: CLIENT_PAGE, offset: intParam('offset', 0),
            state: param('state')};
}

function groupQuery() {
    return {group: state.id,
            group_by_key: param('group_by', 'cluster'),
            state: param('state'),
            limit: CLIENT_PAGE, offset: intParam('offset', 0)};
}

function workunitQuery() {
    return {status: param('status'), task: param('task'),
            assigned_to: param('client'),
            limit: intParam('limit', 50),
            offset: intParam('offset', 0)};
}

function pollInterval() {
    /* Never ask again sooner than a few times what the last round
     * actually cost. On an idle server this is the nominal interval; on
     * a busy one the page quietly gets out of the way. */
    const nominal = state.view === 'overview' ? FAST : SLOW;
    return Math.min(MAX_INTERVAL,
                    Math.max(nominal,
                             BACKOFF_FACTOR * api.lastRoundTripMs));
}

function schedule() {
    if (timer) clearTimeout(timer);
    timer = setTimeout(async () => {
        /* A hidden tab has nobody looking at it; do not poll it. */
        if (!document.hidden) await refresh();
        schedule();
    }, pollInterval());
}

const PAGES = ['overview', 'clients', 'workunits', 'log', 'group',
               'stage'];

function route() {
    const parsed = parseHash();
    const previous = state.view + '/' + state.id;
    state.view = PAGES.includes(parsed.view) ? parsed.view : 'overview';
    state.id = parsed.view === state.view ? parsed.id : null;
    /* A group or a stage is a page about one thing. Named nothing, it
     * is the list it came from, not an empty page of its own. */
    if (state.id === null && TAB_OF[state.view]) {
        state.view = TAB_OF[state.view];
    }
    state.params = parsed.params;
    /* Showing the previous page's detail while the new one loads would
     * be worse than showing nothing. */
    if (previous !== state.view + '/' + state.id) state.detail = null;
    render();
    refresh(true);
}

/* What sort of thing the current page is about, if it is about one
 * thing in particular. */
function detailKind() {
    if (state.id === null) return null;
    if (state.view === 'clients') return 'client';
    if (state.view === 'workunits') return 'workunit';
    if (state.view === 'group') return 'group';
    if (state.view === 'stage') return 'stage';
    return null;
}

/* ---------------- start ---------------- */

async function start() {
    api.readTokenFromFragment();
    const root = document.getElementById('root');

    if (!api.getToken()) {
        clear(root).appendChild(h('main', {},
            tokenGate(async (token) => {
                api.setToken(token);
                try {
                    await api.probe();
                } catch (e) {
                    api.forgetToken();
                    throw e;
                }
                await start();
            })));
        return;
    }

    try {
        await api.probe();
    } catch (e) {
        if (e.status === 401) {
            api.forgetToken();
            return start();
        }
        clear(root).appendChild(h('main', {},
            card('cado-nfs dashboard', banner('error', e.message))));
        return;
    }

    clear(root);
    append(root, [
        h('header', {class: 'topbar', id: 'topbar'}),
        h('nav', {class: 'tabs', id: 'tabs'}),
        h('main', {id: 'main'}),
    ]);
    window.addEventListener('hashchange', route);
    document.addEventListener('visibilitychange', () => {
        if (!document.hidden) refresh(true);
    });
    route();
}

start();
