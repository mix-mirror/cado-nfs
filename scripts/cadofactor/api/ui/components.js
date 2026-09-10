/*
 * DOM helpers and the pieces of the dashboard that more than one view
 * needs. No framework: h() below is the whole of it.
 */

/* Create an element. Children may be nodes, strings, arrays, or null.
 * Attributes starting with "on" become listeners, everything else is
 * set as an attribute -- except for the handful of properties where
 * that would not work. */
export function h(tag, attrs = {}, ...children) {
    const el = document.createElement(tag);
    for (const [key, value] of Object.entries(attrs || {})) {
        if (value === null || value === undefined || value === false) {
            continue;
        }
        if (key.startsWith('on') && typeof value === 'function') {
            el.addEventListener(key.slice(2).toLowerCase(), value);
        } else if (key === 'class') {
            el.className = value;
        } else if (key === 'text') {
            el.textContent = value;
        } else if (key === 'value' || key === 'checked'
                   || key === 'disabled' || key === 'selected') {
            el[key] = value;
        } else {
            el.setAttribute(key, value);
        }
    }
    append(el, children);
    return el;
}

export function append(parent, children) {
    for (const child of children.flat(Infinity)) {
        if (child === null || child === undefined || child === false) {
            continue;
        }
        parent.appendChild(child instanceof Node
                           ? child
                           : document.createTextNode(String(child)));
    }
    return parent;
}

export function clear(el) {
    while (el.firstChild) el.removeChild(el.firstChild);
    return el;
}

/* ---------------- formatting ---------------- */

export function duration(seconds) {
    if (seconds === null || seconds === undefined || isNaN(seconds)) {
        return '–';
    }
    seconds = Math.max(0, seconds);
    /* Keep a decimal below ten seconds: on a small computation a
     * workunit comes back in a fraction of a second, and rounding that
     * to "0s" hides the very number the reader is after. */
    if (seconds < 10) {
        return seconds.toFixed(1).replace(/\.0$/, '') + 's';
    }
    seconds = Math.round(seconds);
    if (seconds < 60) return seconds + 's';
    if (seconds < 3600) {
        return Math.floor(seconds / 60) + 'm ' + (seconds % 60) + 's';
    }
    if (seconds < 86400) {
        return Math.floor(seconds / 3600) + 'h '
            + Math.floor((seconds % 3600) / 60) + 'm';
    }
    return Math.floor(seconds / 86400) + 'd '
        + Math.floor((seconds % 86400) / 3600) + 'h';
}

export function count(value) {
    if (value === null || value === undefined || isNaN(value)) {
        return '–';
    }
    const n = Number(value);
    for (const [limit, suffix] of [[1e9, 'G'], [1e6, 'M'], [1e3, 'k']]) {
        if (Math.abs(n) >= limit) {
            return (n / limit).toPrecision(3) + suffix;
        }
    }
    return String(Math.round(n));
}

export function percent(fraction, digits = 1) {
    if (fraction === null || fraction === undefined || isNaN(fraction)) {
        return '–';
    }
    return (100 * fraction).toFixed(digits) + '%';
}

export function ago(stamp, nowSeconds) {
    if (!stamp) return '–';
    return duration(nowSeconds - stamp);
}

/* ---------------- small pieces ---------------- */

export function pill(text, kind) {
    return h('span', {class: 'pill ' + (kind || 'unknown'), text});
}

export function figure(value, label) {
    return h('div', {class: 'figure'},
             h('div', {class: 'value'}, value),
             h('div', {class: 'label'}, label));
}

export function card(title, ...children) {
    return h('section', {class: 'card'},
             title ? h('h2', {}, title) : null,
             ...children);
}

export function banner(kind, ...children) {
    return h('div', {class: 'banner ' + kind}, ...children);
}

export function empty(message) {
    return h('div', {class: 'empty'}, message);
}

/* A table whose header cells can sort it. `columns` entries are
 * {key, label, num, render, sort}. */
export function table(columns, rows, options = {}) {
    const state = options.state || {};
    const head = h('tr');
    for (const column of columns) {
        const active = state.sort === column.key;
        const arrow = active
            ? h('span', {class: 'arrow'}, state.desc ? ' ▾' : ' ▴')
            : null;
        head.appendChild(h('th', {
            class: (column.num ? 'num ' : '')
                + (column.sort === false ? '' : 'sortable'),
            onclick: column.sort === false ? null : () => {
                if (state.sort === column.key) {
                    state.desc = !state.desc;
                } else {
                    state.sort = column.key;
                    state.desc = !!column.num;
                }
                options.onsort && options.onsort();
            },
        }, column.label, arrow));
    }

    let ordered = rows;
    if (state.sort) {
        const column = columns.find((c) => c.key === state.sort);
        if (column) {
            const key = column.sortValue || ((row) => row[column.key]);
            ordered = rows.slice().sort((a, b) => {
                const x = key(a), y = key(b);
                if (x === y) return 0;
                if (x === null || x === undefined) return 1;
                if (y === null || y === undefined) return -1;
                return (x > y ? 1 : -1) * (state.desc ? -1 : 1);
            });
        }
    }

    const body = h('tbody');
    for (const row of ordered) {
        const tr = h('tr');
        for (const column of columns) {
            const cell = column.render
                ? column.render(row)
                : String(row[column.key] ?? '');
            tr.appendChild(h('td', {
                class: (column.num ? 'num ' : '') + (column.mono ? 'mono' : ''),
            }, cell));
        }
        body.appendChild(tr);
    }

    return h('div', {class: 'tablewrap'},
             h('table', {}, h('thead', {}, head), body));
}

/* ---------------- the token gate ---------------- */

export function tokenGate(onSubmit) {
    const input = h('input', {
        type: 'text',
        placeholder: 'paste the api token here',
        autocomplete: 'off',
        spellcheck: 'false',
    });
    const message = h('div', {});
    const submit = async () => {
        const value = input.value.trim();
        if (!value) return;
        clear(message);
        try {
            await onSubmit(value);
        } catch (e) {
            clear(message).appendChild(banner('error', e.message));
        }
    };
    input.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') submit();
    });
    return h('div', {class: 'card gate'},
             h('h2', {}, 'cado-nfs dashboard'),
             h('p', {},
               'This page needs the api token that the server wrote to ',
               h('code', {}, '<workdir>/<name>.api-token'),
               '. cado-nfs.py also logs a link that carries it, of the ',
               'form ', h('code', {}, '/ui/#token=…'), '.'),
             h('div', {class: 'row'}, input,
               h('button', {class: 'primary', onclick: submit}, 'Open')),
             message);
}
