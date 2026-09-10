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
import {dial, stackedBar, legend, barRows, miniBar} from './charts.js';

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

const VIEWS = [
    {id: 'overview', label: 'Overview'},
    {id: 'clients', label: 'Clients'},
    {id: 'workunits', label: 'Workunits'},
    {id: 'log', label: 'Log'},
];

const state = {
    view: 'overview',
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
    clientSort: {sort: 'completed', desc: true},
    clientPage: {limit: 100, offset: 0},
    wuFilters: {status: '', assigned_to: '', task: '', limit: 50},
    logTail: 200,
};

let timer = null;

/* ---------------- chrome ---------------- */

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
        h('span', {id: 'heartbeat', class: 'heartbeat'}),
        h('button', {
            class: 'small',
            title: 'forget the api token in this tab',
            onclick: () => { api.forgetToken(); location.reload(); },
        }, 'Lock'),
    ]);
    heartbeat();
}

function renderTabs() {
    const nav = document.getElementById('tabs');
    clear(nav);
    for (const view of VIEWS) {
        nav.appendChild(h('a', {
            href: '#' + view.id,
            class: view.id === state.view ? 'active' : '',
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
            count(task.wu_received || 0) + ' / ' + count(task.wu_submitted),
            'workunits back'));
    }
    if (task && task.wu_failed) {
        bits.push(figure(count(task.wu_failed), 'failed'));
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
                    h('div', {class: 'title'}, task.title || task.name),
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
        }
        list.appendChild(h('li', {class: task.phase},
                           h('span', {class: 'dot'}),
                           h('span', {class: 'label'},
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
    }));
}

function workunitsCard() {
    const segments = workunitSegments();
    const total = (state.summary && state.summary.total) || 0;
    return card('Workunits',
                h('div', {class: 'figures'},
                  figure(count(total), 'total'),
                  figure(count((state.summary || {}).outstanding || 0),
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
        .map((k) => h('span', {style: 'margin-right:8px'},
                      pill(counts[k] + ' ' + k, k)));

    const top = rows.map((c) => ({
        label: c.clientid,
        value: c.completed,
        colour: c.state === 'gone' ? '--bad'
            : c.state === 'stale' ? '--warn' : '--accent',
    }));

    return card('Clients',
                h('div', {}, chips, ' ',
                  h('span', {class: 'faint'}, payload.total + ' in all')),
                h('div', {style: 'margin-top:12px'},
                  barRows(top, {format: count})),
                payload.truncated
                    ? h('div', {class: 'faint',
                                style: 'margin-top:8px;font-size:12px'},
                        'and ' + payload.truncated + ' more \u2014 ',
                        h('a', {href: '#clients'}, 'see all'))
                    : null);
}

function strandedCard() {
    const summary = state.clientsSummary || {};
    const counts = summary.counts || {};
    const quiet = (counts.stale || 0) + (counts.gone || 0);
    if (!quiet) return null;
    /* Naming them needs the full list, which the overview does not
     * poll. Point at the clients view instead of dragging it in. */
    const rows = ((state.clients || {}).clients || [])
        .filter((c) => c.in_flight
                && (c.state === 'stale' || c.state === 'gone'));
    if (!rows.length) {
        return h('section', {class: 'card span2'},
                 banner('warn',
                        h('strong', {}, quiet + ' client'
                          + (quiet === 1 ? '' : 's')),
                        ' have gone quiet. ',
                        h('a', {href: '#clients'},
                          'Look at them on the clients page'),
                        ' to reclaim what they are holding.'));
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
                 {key: 'clientid', label: 'client', mono: true},
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
                   onclick: (e) => bulkReclaim(e.target),
               }, 'Reclaim everything assigned longer ago than'),
               h('input', {type: 'number', id: 'bulk-minutes', value: 60,
                           min: 1, style: 'width:80px'}),
               h('span', {class: 'muted'}, 'minutes')));
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
    ]);
}

/* ---------------- clients ---------------- */

function clientsView() {
    const main = document.getElementById('main');
    clear(main);
    if (state.notice) main.appendChild(state.notice);
    const payload = state.clients || {};
    const rows = payload.clients || [];
    if (!rows.length) {
        main.appendChild(card('Clients',
                              empty('No client has asked for work yet.')));
        return;
    }
    const totalDone = rows.reduce((s, c) => s + c.completed, 0) || 1;
    const page = state.clientPage;
    const total = payload.total || rows.length;
    const pager = total > rows.length
        ? h('div', {class: 'controls', style: 'margin-bottom:0'},
            h('button', {
                class: 'small', disabled: page.offset <= 0,
                onclick: () => {
                    page.offset = Math.max(0, page.offset - page.limit);
                    refresh(true);
                },
            }, '\u2190 previous'),
            h('span', {class: 'muted'},
              (page.offset + 1) + '\u2013'
              + (page.offset + rows.length) + ' of ' + total),
            h('button', {
                class: 'small',
                disabled: page.offset + rows.length >= total,
                onclick: () => {
                    page.offset += page.limit;
                    refresh(true);
                },
            }, 'next \u2192'))
        : null;
    main.appendChild(card(
        'Clients',
        h('p', {class: 'muted', style: 'margin-top:-6px'},
          'How long a client may be silent before it counts as stale is ',
          'judged per client, from how long its own workunits have ',
          'recently been taking \u2014 hover a state to see the ',
          'reasoning. It is never longer than tasks.wutimeout (',
          duration(payload.wutimeout), '), which is when its work gets ',
          'reassigned anyway, and while few workunits back the ',
          'estimate it is held down to a couple of ',
          'tasks.wutimeoutcheck intervals.'),
        pager,
        table([
            {key: 'clientid', label: 'client', mono: true},
            {key: 'state', label: 'state',
             render: (r) => statePill(r)},
            {key: 'typical_turnaround', label: 'usual pace', num: true,
             render: (r) => r.typical_turnaround === null
                 ? h('span', {class: 'faint',
                              title: 'not enough recent workunits'
                                     + ' from this client yet'}, '–')
                 : h('span', {title: 'median of '
                                     + r.turnaround_samples
                                     + ' recent workunits'},
                     duration(r.typical_turnaround))},
            {key: 'in_flight', label: 'in flight', num: true},
            {key: 'completed', label: 'completed', num: true},
            {key: 'failed', label: 'failed', num: true,
             render: (r) => r.failed
                 ? h('span', {style: 'color:var(--bad)'}, r.failed)
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
        ], rows, {state: state.clientSort, onsort: () => clientsView()}),
        pager ? h('div', {class: 'faint',
                          style: 'margin-top:10px;font-size:12px'},
                  'Sorting applies to this page only: the server orders'
                  + ' clients by what they have contributed, and the'
                  + ' dashboard asks for one page at a time so that a'
                  + ' large pool does not cost the computation'
                  + ' bandwidth it needs for workunits.')
              : null));
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

    const filters = h('div', {class: 'controls'},
        h('label', {}, 'status',
          h('select', {
              onchange: (e) => {
                  state.wuFilters.status = e.target.value;
                  refresh(true);
              },
          }, STATUSES.map((s) => h('option', {
              value: s, selected: s === state.wuFilters.status,
          }, s || 'any')))),
        h('label', {}, 'task',
          h('select', {
              onchange: (e) => {
                  state.wuFilters.task = e.target.value;
                  refresh(true);
              },
          }, [''].concat(tasks).map((t) => h('option', {
              value: t, selected: t === state.wuFilters.task,
          }, t || 'any')))),
        h('label', {}, 'client',
          h('input', {
              type: 'text', value: state.wuFilters.assigned_to,
              placeholder: 'assigned to',
              onchange: (e) => {
                  state.wuFilters.assigned_to = e.target.value.trim();
                  refresh(true);
              },
          })),
        h('label', {}, 'rows',
          h('select', {
              onchange: (e) => {
                  state.wuFilters.limit = Number(e.target.value);
                  refresh(true);
              },
          }, [25, 50, 100, 200].map((n) => h('option', {
              value: n, selected: n === state.wuFilters.limit,
          }, n)))));

    main.appendChild(card('Workunits', filters,
        rows.length
            ? table([
                {key: 'wuid', label: 'workunit', mono: true},
                {key: 'status_name', label: 'status',
                 render: (r) => statusPill(r.status_name)},
                {key: 'task', label: 'task'},
                {key: 'attempt', label: 'try', num: true},
                {key: 'assignedclient', label: 'client',
                 render: (r) => r.assignedclient || r.resultclient
                     || h('span', {class: 'faint'}, '–')},
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
                onchange: (e) => {
                    state.logTail = Number(e.target.value);
                    refresh(true);
                },
            }, [100, 200, 500, 1000, 2000].map((n) => h('option', {
                value: n, selected: n === state.logTail,
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
    ({
        overview,
        clients: clientsView,
        workunits: workunitsView,
        log: logView,
    })[state.view]();
}

async function refresh(immediate = false) {
    const started = Date.now();
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
        if (state.view === 'overview') {
            wanted.push(api.summary().then((r) => { state.summary = r; }));
            wanted.push(api.clientsSummary()
                .then((r) => { state.clientsSummary = r; }));
        }
        if (state.view === 'clients') {
            wanted.push(api.clientsSummary()
                .then((r) => { state.clientsSummary = r; }));
            wanted.push(api.clients(state.clientPage)
                .then((r) => { state.clients = r; }));
        }
        if (state.view === 'workunits') {
            wanted.push(api.summary().then((r) => { state.summary = r; }));
            wanted.push(api.workunits(state.wuFilters)
                .then((r) => { state.workunits = r; }));
        }
        if (state.view === 'log') {
            wanted.push(api.log(state.logTail)
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
        state.error = e.message;
    }
    api.noteRoundTrip(Date.now() - started);
    render();
    if (immediate) schedule();
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

function route() {
    const id = (location.hash || '#overview').slice(1);
    state.view = VIEWS.some((v) => v.id === id) ? id : 'overview';
    render();
    refresh(true);
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
