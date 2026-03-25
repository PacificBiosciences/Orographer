import contextlib
import http.client
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest


def _get_free_port() -> int:
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.mark.slow
def test_deploy_smoke_serves_html():
    # Use an existing expected_outputs directory so deploy has static content.
    outdir = Path("tests/data/expected_outputs/paraviewer/HG002.paraphase.small/hba").resolve()
    assert outdir.is_dir()

    port = _get_free_port()
    server_cmd = [
        sys.executable,
        "-m",
        "orographer",
        "deploy",
        "--outdir",
        str(outdir),
        "--port",
        str(port),
    ]

    # Start server
    proc = subprocess.Popen(
        server_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        deadline = time.time() + 30
        url_path = "/chr16_171801_175500_bokeh.html"
        while time.time() < deadline:
            # If the server exited early, surface stderr to help debugging.
            if proc.poll() is not None:
                stderr = ""
                with contextlib.suppress(Exception):
                    stderr = (proc.stderr.read() or "").strip()  # type: ignore[union-attr]
                raise AssertionError(f"Deploy server exited early. stderr={stderr!r}")

            try:
                conn = http.client.HTTPConnection("127.0.0.1", port, timeout=1)
                conn.request("GET", url_path)
                resp = conn.getresponse()
                body = resp.read()
                if resp.status != 200:
                    raise AssertionError(f"Expected 200, got {resp.status}")
                assert b"Bokeh" in body or b"root\\u002EBokeh" in body
                return
            except (ConnectionRefusedError, TimeoutError):
                time.sleep(0.25)

        raise AssertionError("Deploy server did not respond in time")
    finally:
        proc.terminate()
        with contextlib.suppress(Exception):
            proc.wait(timeout=10)
        with contextlib.suppress(Exception):
            if proc.poll() is None:
                proc.kill()
