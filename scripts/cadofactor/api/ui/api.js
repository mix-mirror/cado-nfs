/*
 * Talking to the cado-nfs api from the browser.
 *
 * Everything the dashboard shows comes through here, and every call
 * below corresponds to an operation in /api/v1/openapi.json. There is
 * deliberately no other path to the data.
 */

const PREFIX = '/api/v1';

/* Where the token lives for the length of this tab. sessionStorage
 * rather than localStorage: a token is a credential, and it has no
 * business outliving the window it was pasted into. */
const TOKEN_KEY = 'cado-api-token';

export class ApiError extends Error {
    constructor(status, message) {
        super(message);
        this.status = status;
    }
}

/* Kept out of the response bodies on purpose -- see json_response() in
 * admin.py -- so that unchanged answers can revalidate as 304. */
let serverTime = null;
let serverTimeSeenAt = null;

/* Cached bodies, keyed by url, so a 304 costs nothing to handle. */
const etags = new Map();
const bodies = new Map();

export function readTokenFromFragment() {
    /* cado-nfs.py prints a url ending in #token=... . A fragment is
     * never sent to the server, so it appears in no access log and no
     * Referer header. Move it into session storage and wipe it from the
     * address bar so it does not linger in history or in a screenshot. */
    const match = /(?:^|[#&])token=([^&]+)/.exec(window.location.hash || '');
    if (!match) return null;
    const token = decodeURIComponent(match[1]);
    setToken(token);
    history.replaceState(null, '', window.location.pathname +
                         window.location.search);
    return token;
}

export function getToken() {
    try {
        return window.sessionStorage.getItem(TOKEN_KEY);
    } catch (e) {
        /* Private mode, or storage blocked. Hold it in memory instead. */
        return memoryToken;
    }
}

let memoryToken = null;

export function setToken(token) {
    memoryToken = token;
    try {
        window.sessionStorage.setItem(TOKEN_KEY, token);
    } catch (e) { /* memoryToken is the fallback */ }
}

export function forgetToken() {
    memoryToken = null;
    try {
        window.sessionStorage.removeItem(TOKEN_KEY);
    } catch (e) { /* nothing to do */ }
    etags.clear();
    bodies.clear();
}

/* The server's clock, advanced by however long ago we last heard it.
 * Ages are computed from this rather than from Date.now(), so that a
 * skewed browser clock cannot make a healthy client look stale. */
export function now() {
    if (serverTime === null) return Date.now() / 1000;
    return serverTime + (Date.now() - serverTimeSeenAt) / 1000;
}

function noteServerTime(response) {
    const stamp = response.headers.get('X-Cado-Server-Time');
    if (stamp) {
        const value = parseFloat(stamp);
        if (!isNaN(value)) {
            serverTime = value;
            serverTimeSeenAt = Date.now();
        }
    }
}

async function describe(response) {
    try {
        const body = await response.json();
        if (body && body.description) return body.description;
    } catch (e) { /* not JSON, fall through */ }
    return response.statusText || ('HTTP ' + response.status);
}

async function call(path, {method = 'GET', body = null,
                           revalidate = false} = {}) {
    const headers = {'Accept': 'application/json'};
    const token = getToken();
    if (token) headers['Authorization'] = 'Bearer ' + token;
    if (body !== null) headers['Content-Type'] = 'application/json';
    if (revalidate && etags.has(path)) {
        headers['If-None-Match'] = etags.get(path);
    }

    let response;
    try {
        response = await fetch(path, {
            method,
            headers,
            body: body === null ? undefined : JSON.stringify(body),
            cache: 'no-store',
        });
    } catch (e) {
        throw new ApiError(0, 'cannot reach the server (' + e.message + ')');
    }
    noteServerTime(response);

    if (response.status === 304 && bodies.has(path)) {
        return bodies.get(path);
    }
    if (!response.ok) {
        throw new ApiError(response.status, await describe(response));
    }

    const etag = response.headers.get('ETag');
    const payload = response.status === 204 ? null : await response.json();
    if (etag && method === 'GET') {
        etags.set(path, etag);
        bodies.set(path, payload);
    }
    return payload;
}

/* Reads. Each of these polls, so each revalidates. */
export const info = () => call(PREFIX + '/info', {revalidate: true});
export const progress = () => call(PREFIX + '/progress', {revalidate: true});
export const summary = () =>
    call(PREFIX + '/workunits/summary', {revalidate: true});
/* The tallies and the top contributors only. This is what the overview
 * polls: the full list is half a megabyte on a 1400-client run, and
 * asking for that every two seconds would take capacity away from the
 * computation being watched. */
export const clientsSummary = () =>
    call(PREFIX + '/clients?summary=1', {revalidate: true});

export function clients(query = {}) {
    const params = new URLSearchParams();
    for (const [key, value] of Object.entries(query)) {
        if (value !== null && value !== undefined && value !== '') {
            params.set(key, value);
        }
    }
    const q = params.toString();
    return call(PREFIX + '/clients' + (q ? '?' + q : ''),
                {revalidate: true});
}

export const stats = () => call(PREFIX + '/stats', {revalidate: true});
export const parameters = () =>
    call(PREFIX + '/parameters', {revalidate: true});

export async function setParameter(name, value) {
    const r = await call(PREFIX + '/parameters/' + encodeURIComponent(name),
                         {method: 'POST', body: {value}});
    invalidate();
    return r;
}

export function workunits(query = {}) {
    const params = new URLSearchParams();
    for (const [key, value] of Object.entries(query)) {
        if (value !== null && value !== undefined && value !== '') {
            params.set(key, value);
        }
    }
    const q = params.toString();
    return call(PREFIX + '/workunits' + (q ? '?' + q : ''),
                {revalidate: true});
}

export const workunit = (wuid) =>
    call(PREFIX + '/workunits/' + encodeURIComponent(wuid),
         {revalidate: true});

export const client = (clientid) =>
    call(PREFIX + '/clients/' + encodeURIComponent(clientid),
         {revalidate: true});

export const log = (tail) =>
    call(PREFIX + '/log?tail=' + encodeURIComponent(tail),
         {revalidate: true});

/* Writes. Never revalidated, and they invalidate what they change. */

function invalidate() {
    etags.clear();
    bodies.clear();
}

export async function reclaimClient(clientid) {
    const r = await call(
        PREFIX + '/clients/' + encodeURIComponent(clientid) + '/reclaim',
        {method: 'POST'});
    invalidate();
    return r;
}

export async function resubmitWorkunit(wuid) {
    const r = await call(
        PREFIX + '/workunits/' + encodeURIComponent(wuid) + '/resubmit',
        {method: 'POST'});
    invalidate();
    return r;
}

/* The group forms. One request per group rather than one per client:
 * a cluster can be a few hundred of them. */
export async function reclaimClients(selector) {
    const r = await call(PREFIX + '/clients/reclaim',
                         {method: 'POST', body: selector});
    invalidate();
    return r;
}

export async function reclaimOlderThan(seconds) {
    const r = await call(PREFIX + '/workunits/reclaim',
                         {method: 'POST', body: {older_than: seconds}});
    invalidate();
    return r;
}

export async function setServing(serving) {
    const r = await call(PREFIX + '/serving',
                         {method: 'POST', body: {serving}});
    invalidate();
    return r;
}

/* Used by the token gate to find out whether a token works at all. */
export async function probe() {
    await call(PREFIX + '/info');
    return true;
}

/* How long the last round of polling took. The dashboard uses it to
 * back off: a server that is slow to answer is a server busy handing
 * out workunits, and a monitor has no business elbowing in. */
export let lastRoundTripMs = 0;

export function noteRoundTrip(ms) {
    lastRoundTripMs = ms;
}
