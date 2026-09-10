/*
 * Charts, drawn as inline SVG by hand.
 *
 * There is no charting library here and there is not meant to be one:
 * the three shapes this dashboard needs are a progress dial, a stacked
 * bar and a row of horizontal bars, and all three are a few lines of
 * geometry. Colours come from the CSS custom properties in style.css,
 * so light and dark mode need no work here.
 */

const NS = 'http://www.w3.org/2000/svg';

function svgel(tag, attrs = {}, ...children) {
    const el = document.createElementNS(NS, tag);
    for (const [key, value] of Object.entries(attrs)) {
        if (value === null || value === undefined) continue;
        el.setAttribute(key, String(value));
    }
    for (const child of children.flat(Infinity)) {
        if (child === null || child === undefined) continue;
        el.appendChild(child instanceof Node
                       ? child
                       : document.createTextNode(String(child)));
    }
    return el;
}

function css(name) {
    return getComputedStyle(document.documentElement)
        .getPropertyValue(name).trim() || '#888';
}

/* A ring showing one fraction, with a caption in the middle. Used for
 * the achievement of the task that is running. */
export function dial(fraction, caption, subcaption, size = 128) {
    const stroke = 11;
    const r = (size - stroke) / 2;
    const c = size / 2;
    const circumference = 2 * Math.PI * r;
    const value = fraction === null || fraction === undefined
        ? 0 : Math.max(0, Math.min(1, fraction));

    const svg = svgel('svg', {
        class: 'chart',
        width: size,
        height: size,
        viewBox: `0 0 ${size} ${size}`,
        role: 'img',
        'aria-label': caption + ' complete',
    },
        svgel('circle', {
            cx: c, cy: c, r,
            fill: 'none',
            stroke: css('--border'),
            'stroke-width': stroke,
        }),
        svgel('circle', {
            cx: c, cy: c, r,
            fill: 'none',
            stroke: css('--accent'),
            'stroke-width': stroke,
            'stroke-linecap': 'round',
            'stroke-dasharray': `${circumference * value} ${circumference}`,
            /* Start at twelve o'clock rather than at three. */
            transform: `rotate(-90 ${c} ${c})`,
        }),
        svgel('text', {
            x: c, y: c - 2,
            'text-anchor': 'middle',
            'dominant-baseline': 'middle',
            fill: css('--text'),
            'font-size': 21,
            'font-weight': 600,
        }, caption),
        subcaption ? svgel('text', {
            x: c, y: c + 18,
            'text-anchor': 'middle',
            'dominant-baseline': 'middle',
            fill: css('--faint'),
            'font-size': 10.5,
        }, subcaption) : null);
    return svg;
}

/* One bar split into labelled segments. `segments` are
 * {label, value, colour} where colour names a CSS custom property. */
export function stackedBar(segments, width = 520, height = 26) {
    const total = segments.reduce((sum, s) => sum + s.value, 0);
    const svg = svgel('svg', {
        class: 'chart',
        viewBox: `0 0 ${width} ${height}`,
        preserveAspectRatio: 'none',
        height,
        role: 'img',
        'aria-label': segments.map((s) => `${s.value} ${s.label}`)
            .join(', '),
    });
    if (!total) {
        svg.appendChild(svgel('rect', {
            x: 0, y: 0, width, height, rx: 5,
            fill: css('--border'),
        }));
        return svg;
    }
    let x = 0;
    for (const segment of segments) {
        if (!segment.value) continue;
        const w = (segment.value / total) * width;
        svg.appendChild(svgel('rect', {
            x, y: 0, width: Math.max(w, 1.5), height,
            fill: css(segment.colour),
        }, svgel('title', {},
                 `${segment.label}: ${segment.value}`)));
        x += w;
    }
    /* Rounded ends without clipping the segment colours. */
    svg.setAttribute('style', 'border-radius:5px;overflow:hidden');
    return svg;
}

export function legend(segments) {
    const wrap = document.createElement('div');
    wrap.className = 'legend';
    for (const segment of segments) {
        if (segment.value === 0 && segment.hideEmpty) continue;
        const item = document.createElement('span');
        const swatch = document.createElement('span');
        swatch.className = 'swatch';
        swatch.style.background = css(segment.colour);
        item.appendChild(swatch);
        item.appendChild(document.createTextNode(
            `${segment.label} ${segment.value}`));
        wrap.appendChild(item);
    }
    return wrap;
}

/* Horizontal bars, one per row: {label, value, colour}. Used for what
 * each client has contributed. */
export function barRows(rows, {width = 520, rowHeight = 22,
                               labelWidth = 130,
                               format = (v) => String(v)} = {}) {
    const max = rows.reduce((m, r) => Math.max(m, r.value), 0) || 1;
    const height = Math.max(rows.length * rowHeight, rowHeight);
    const barWidth = width - labelWidth - 62;
    const svg = svgel('svg', {
        class: 'chart',
        viewBox: `0 0 ${width} ${height}`,
        height,
        role: 'img',
        'aria-label': 'contribution per client',
    });
    rows.forEach((row, i) => {
        const y = i * rowHeight;
        const w = Math.max((row.value / max) * barWidth, row.value ? 2 : 0);
        svg.appendChild(svgel('text', {
            x: 0, y: y + rowHeight / 2,
            'dominant-baseline': 'middle',
            fill: css('--text'),
            'font-size': 12,
        }, clip(row.label, 18)));
        svg.appendChild(svgel('rect', {
            x: labelWidth, y: y + 4,
            width: w, height: rowHeight - 9,
            rx: 3,
            fill: css(row.colour || '--accent'),
        }, svgel('title', {}, `${row.label}: ${format(row.value)}`)));
        svg.appendChild(svgel('text', {
            x: labelWidth + w + 7, y: y + rowHeight / 2,
            'dominant-baseline': 'middle',
            fill: css('--muted'),
            'font-size': 11.5,
        }, format(row.value)));
    });
    return svg;
}

function clip(text, n) {
    return text.length > n ? text.slice(0, n - 1) + '…' : text;
}

/* A slim bar with no labels, for use inside a table cell. */
export function miniBar(fraction, width = 70, height = 7) {
    const value = Math.max(0, Math.min(1, fraction || 0));
    /* No 'chart' class here: that makes it display:block, which would
     * push whatever sits next to it in a table cell onto its own line. */
    return svgel('svg', {
        width, height,
        viewBox: `0 0 ${width} ${height}`,
        style: 'vertical-align:middle',
    },
        svgel('rect', {
            x: 0, y: 0, width, height, rx: height / 2,
            fill: css('--border'),
        }),
        svgel('rect', {
            x: 0, y: 0, width: value * width, height, rx: height / 2,
            fill: css('--accent'),
        }));
}
