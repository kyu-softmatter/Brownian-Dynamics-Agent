"""Minimal line-coverage for bdbot/, via sys.monitoring. A pytest plugin."""
import sys, os, json, atexit
MON = sys.monitoring
TOOL = 2
ROOT = os.path.abspath("bdbot")
hits = {}
def _line(code, lineno):
    fn = code.co_filename
    if fn.startswith(ROOT):
        hits.setdefault(fn, set()).add(lineno)
    return MON.DISABLE if not fn.startswith(ROOT) else None
def pytest_configure(config):
    MON.use_tool_id(TOOL, "bdcov")
    MON.register_callback(TOOL, MON.events.LINE, _line)
    MON.set_events(TOOL, MON.events.LINE)
def pytest_unconfigure(config):
    try:
        MON.set_events(TOOL, 0); MON.free_tool_id(TOOL)
    except Exception: pass
    json.dump({k: sorted(v) for k, v in hits.items()}, open("/tmp/bdcov.json","w"))
