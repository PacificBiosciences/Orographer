from __future__ import annotations

from pathlib import Path

import pytest

from orographer import deploy


def test_run_deploy_exits_when_output_dir_is_missing(tmp_path: Path) -> None:
    missing_dir = tmp_path / "missing"

    with pytest.raises(SystemExit) as exc_info:
        deploy.run_deploy(str(missing_dir), port=8765)

    assert exc_info.value.code == 1


def test_run_deploy_changes_to_output_dir_and_serves(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    calls: dict[str, object] = {}

    class FakeTCPServer:
        def __init__(self, server_address: tuple[str, int], handler_class: type) -> None:
            calls["server_address"] = server_address
            calls["handler_class"] = handler_class

        def __enter__(self) -> FakeTCPServer:
            return self

        def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
            calls["exited"] = True

        def serve_forever(self) -> None:
            calls["served"] = True

    monkeypatch.setattr(deploy.os, "chdir", lambda path: calls.setdefault("chdir", path))
    monkeypatch.setattr(deploy.socketserver, "TCPServer", FakeTCPServer)

    deploy.run_deploy(str(tmp_path), port=8765)

    captured = capsys.readouterr()
    assert calls["chdir"] == str(tmp_path)
    assert calls["server_address"] == ("", 8765)
    assert calls["handler_class"] is deploy.PlotHTTPRequestHandler
    assert calls["served"] is True
    assert calls["exited"] is True
    assert "Serving plots from:" in captured.out


def test_plot_request_handler_serves_static_files_only() -> None:
    assert issubclass(deploy.PlotHTTPRequestHandler, deploy.http.server.SimpleHTTPRequestHandler)
    assert "do_GET" not in deploy.PlotHTTPRequestHandler.__dict__
