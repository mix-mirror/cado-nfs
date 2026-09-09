"""
OpenAPI 3.1 description of the cado-nfs api.

The document is assembled here, at run time, from metadata that each
endpoint carries. No third-party tooling is involved -- no apispec, no
marshmallow, no flasgger -- so `GET /api/v1/openapi.json` answers on any
installation that can run the server at all. That is deliberate: the
previous arrangement embedded Swagger 2.0 fragments in docstrings and
produced no machine-readable description whatsoever when flasgger
happened to be absent.

Routes and documentation come from a single declaration:

    @api_route("/api/v1/clients", tags=["monitoring"], auth=True,
               summary="Clients seen by the server",
               responses={200: ("Client list", CLIENT_LIST)})
    def api_clients(self):
        ...

`collect_api_routes()` finds the decorated methods, the caller registers
them with flask, and `build_spec()` turns the same metadata into the
OpenAPI document. An endpoint cannot therefore be routed but undocumented
-- it would have to be registered by hand to escape, which is what
tests/scripts/cadofactor/test_api_openapi.sh checks for.
"""

import re

OPENAPI_VERSION = "3.1.0"

# Placeholders in a werkzeug rule, e.g. <wuid> or <path:path>.
_RULE_ARG = re.compile(r"<(?:([a-zA-Z_][a-zA-Z0-9_]*):)?"
                       r"([a-zA-Z_][a-zA-Z0-9_]*)>")

_CONVERTER_SCHEMA = {
    "string": {"type": "string"},
    "path": {"type": "string"},
    "int": {"type": "integer"},
    "float": {"type": "number"},
    "uuid": {"type": "string", "format": "uuid"},
}

# Methods that werkzeug adds on its own and that carry no documentation.
_IMPLICIT_METHODS = ("HEAD", "OPTIONS")


def rule_to_openapi_path(rule):
    """
    Translate a werkzeug rule into an OpenAPI path template.

    >>> rule_to_openapi_path('/WUstatus/<wuid>')
    '/WUstatus/{wuid}'
    >>> rule_to_openapi_path('/file/<path:path>')
    '/file/{path}'
    >>> rule_to_openapi_path('/api/v1/clients/<clientid>/reclaim')
    '/api/v1/clients/{clientid}/reclaim'
    >>> rule_to_openapi_path('/')
    '/'
    """
    return _RULE_ARG.sub(lambda m: "{%s}" % m.group(2), rule)


def rule_parameters(rule):
    """
    Build the OpenAPI parameter objects implied by a rule's placeholders.

    >>> [p['name'] for p in rule_parameters('/a/<int:n>/b/<path:p>')]
    ['n', 'p']
    >>> rule_parameters('/a/<int:n>')[0]['schema']['type']
    'integer'
    >>> rule_parameters('/a/<n>')[0]['schema']['type']
    'string'
    >>> rule_parameters('/a/<n>')[0]['required']
    True
    >>> rule_parameters('/plain')
    []
    """
    out = []
    for m in _RULE_ARG.finditer(rule):
        converter, name = m.group(1), m.group(2)
        schema = _CONVERTER_SCHEMA.get(converter, {"type": "string"})
        out.append({"name": name,
                    "in": "path",
                    "required": True,
                    "schema": dict(schema)})
    return out


def query_parameter(name, schema, description, required=False):
    """
    Shorthand for a query string parameter object.

    >>> p = query_parameter('limit', {'type': 'integer'}, 'how many')
    >>> p['in'], p['name'], p['required']
    ('query', 'limit', False)
    """
    return {"name": name,
            "in": "query",
            "required": required,
            "description": description,
            "schema": schema}


def _normalize_responses(responses):
    """
    Accept a terse spelling of the responses map and expand it.

    A value may be a bare description, or a (description, schema) pair.

    >>> r = _normalize_responses({200: ('ok', {'type': 'object'}),
    ...                           404: 'nope'})
    >>> sorted(r)
    ['200', '404']
    >>> r['404']
    {'description': 'nope'}
    >>> r['200']['content']['application/json']['schema']
    {'type': 'object'}
    """
    out = {}
    for code, value in (responses or {}).items():
        if isinstance(value, tuple):
            description, schema = value
        else:
            description, schema = value, None
        entry = {"description": description}
        if schema is not None:
            entry["content"] = {"application/json": {"schema": schema}}
        out[str(code)] = entry
    return out


def json_body(schema, description=None, required=True):
    """
    Shorthand for a JSON request body object.

    >>> b = json_body({'type': 'object'})
    >>> b['required']
    True
    >>> 'application/json' in b['content']
    True
    """
    body = {"required": required,
            "content": {"application/json": {"schema": schema}}}
    if description is not None:
        body["description"] = description
    return body


def api_route(rule, methods=("GET",), *,
              tags=None, summary=None, description=None,
              parameters=None, request_body=None, responses=None,
              auth=False, produces_json=True):
    """
    Declare an endpoint: its route *and* its documentation.

    `auth` marks the endpoint as requiring the bearer token. It only
    records the fact in the document; enforcement is the business of
    cadofactor.api.auth.require_token.

    >>> class T:
    ...     @api_route("/x/<n>", summary="hi", tags=["misc"])
    ...     def api_x(self, n):
    ...         return n
    >>> T.api_x.cado_api_route['rule']
    '/x/<n>'
    >>> T.api_x.cado_api_route['methods']
    ('GET',)
    >>> T().api_x(3)
    3
    """
    def decorate(fn):
        operation = {"operationId": fn.__name__,
                     "tags": list(tags or ["misc"])}
        if summary is not None:
            operation["summary"] = summary
        if description is not None:
            operation["description"] = description
        params = rule_parameters(rule) + list(parameters or [])
        if params:
            operation["parameters"] = params
        if request_body is not None:
            operation["requestBody"] = request_body
        if responses:
            operation["responses"] = _normalize_responses(responses)
        elif produces_json:
            operation["responses"] = _normalize_responses(
                {200: "Success"})
        if auth:
            operation["security"] = [{"bearerAuth": []}]
        fn.cado_api_route = {"rule": rule,
                             "methods": tuple(methods),
                             "auth": auth,
                             "operation": operation}
        return fn
    return decorate


def collect_api_routes(obj):
    """
    Return the api endpoints declared on obj, sorted by rule.

    Each element is a (bound_method, metadata) pair. Attributes are read
    off the *class*, so that properties are not evaluated in passing.

    >>> class T:
    ...     @api_route("/b")
    ...     def api_b(self):
    ...         pass
    ...     @api_route("/a")
    ...     def api_a(self):
    ...         pass
    ...     def helper(self):
    ...         pass
    >>> [m['rule'] for _, m in collect_api_routes(T())]
    ['/a', '/b']
    """
    found = []
    cls = type(obj)
    for name in dir(cls):
        attribute = getattr(cls, name, None)
        meta = getattr(attribute, "cado_api_route", None)
        if meta is not None:
            found.append((meta["rule"], name, getattr(obj, name), meta))
    found.sort(key=lambda t: (t[0], t[1]))
    return [(method, meta) for _, _, method, meta in found]


def build_spec(routes, *, title, version, description=None, servers=None,
               license_name=None, license_url=None, contact=None):
    """
    Assemble the OpenAPI document from collected route metadata.

    >>> class T:
    ...     @api_route("/api/v1/x", methods=("GET", "POST"), auth=True,
    ...                summary="s")
    ...     def api_x(self):
    ...         pass
    >>> doc = build_spec([m for _, m in collect_api_routes(T())],
    ...                  title='t', version='1')
    >>> doc['openapi']
    '3.1.0'
    >>> sorted(doc['paths']['/api/v1/x'])
    ['get', 'post']
    >>> doc['paths']['/api/v1/x']['get']['security']
    [{'bearerAuth': []}]
    >>> doc['components']['securitySchemes']['bearerAuth']['scheme']
    'bearer'
    """
    paths = {}
    used_tags = set()
    for meta in routes:
        path = rule_to_openapi_path(meta["rule"])
        entry = paths.setdefault(path, {})
        used_tags.update(meta["operation"].get("tags", []))
        for method in meta["methods"]:
            if method.upper() in _IMPLICIT_METHODS:
                continue
            operation = dict(meta["operation"])
            # operationId must be unique across the document, and a
            # single view function may serve several methods.
            if len(meta["methods"]) > 1:
                operation["operationId"] = "%s_%s" % (
                    operation["operationId"], method.lower())
            entry[method.lower()] = operation

    info = {"title": title, "version": version}
    if description is not None:
        info["description"] = description
    if license_name is not None:
        info["license"] = {"name": license_name}
        if license_url is not None:
            info["license"]["url"] = license_url
    if contact is not None:
        info["contact"] = contact

    doc = {
        "openapi": OPENAPI_VERSION,
        "info": info,
        "paths": paths,
        "components": {
            "securitySchemes": {
                "bearerAuth": {
                    "type": "http",
                    "scheme": "bearer",
                    "description":
                        "Token written by the server to"
                        " <workdir>/<name>.api-token with mode 0600."
                        " Send it as 'Authorization: Bearer <token>'.",
                },
            },
        },
        "tags": [{"name": t} for t in sorted(used_tags)],
    }
    if servers:
        doc["servers"] = [{"url": u} for u in servers]
    return doc


if __name__ == "__main__":
    import doctest
    doctest.testmod()
