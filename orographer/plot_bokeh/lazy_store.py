"""Reusable static payload writers for lazy-loaded Bokeh plot data."""

import gzip
import json
from typing import Any


def write_json_gzip(path: str, payload: Any) -> None:
    """Write compact JSON in gzip format for static deployment."""
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6) as handle:
        json.dump(payload, handle, separators=(",", ":"))


def _pack_int(value: int) -> bytes:
    if 0 <= value <= 0x7F:
        return bytes([value])
    if 0 <= value <= 0xFF:
        return b"\xcc" + value.to_bytes(1, "big")
    if 0 <= value <= 0xFFFF:
        return b"\xcd" + value.to_bytes(2, "big")
    if 0 <= value <= 0xFFFFFFFF:
        return b"\xce" + value.to_bytes(4, "big")
    if 0 <= value <= 0xFFFFFFFFFFFFFFFF:
        return b"\xcf" + value.to_bytes(8, "big")
    if -32 <= value < 0:
        return bytes([0x100 + value])
    if -128 <= value < 0:
        return b"\xd0" + value.to_bytes(1, "big", signed=True)
    if -32768 <= value < 0:
        return b"\xd1" + value.to_bytes(2, "big", signed=True)
    if -2147483648 <= value < 0:
        return b"\xd2" + value.to_bytes(4, "big", signed=True)
    return b"\xd3" + value.to_bytes(8, "big", signed=True)


def _pack_sequence_header(prefix: int, code16: bytes, code32: bytes, length: int) -> bytes:
    if length <= 15:
        return bytes([prefix + length])
    if length <= 0xFFFF:
        return code16 + length.to_bytes(2, "big")
    return code32 + length.to_bytes(4, "big")


def pack_msgpack(value: Any) -> bytes:
    """Pack the small MessagePack subset used by static lazy-load shards."""
    if value is None:
        return b"\xc0"
    if value is False:
        return b"\xc2"
    if value is True:
        return b"\xc3"
    if isinstance(value, int):
        return _pack_int(value)
    if isinstance(value, str):
        encoded = value.encode("utf-8")
        length = len(encoded)
        if length <= 31:
            return bytes([0xA0 + length]) + encoded
        if length <= 0xFF:
            return b"\xd9" + length.to_bytes(1, "big") + encoded
        if length <= 0xFFFF:
            return b"\xda" + length.to_bytes(2, "big") + encoded
        return b"\xdb" + length.to_bytes(4, "big") + encoded
    if isinstance(value, list | tuple):
        header = _pack_sequence_header(0x90, b"\xdc", b"\xdd", len(value))
        return header + b"".join(pack_msgpack(item) for item in value)
    if isinstance(value, dict):
        header = _pack_sequence_header(0x80, b"\xde", b"\xdf", len(value))
        return header + b"".join(
            pack_msgpack(key) + pack_msgpack(item) for key, item in value.items()
        )
    raise TypeError(f"Unsupported MessagePack value: {type(value)!r}")


def write_msgpack_gzip(path: str, payload: Any) -> None:
    """Write compact MessagePack in gzip format for static deployment."""
    with gzip.open(path, "wb", compresslevel=6) as handle:
        handle.write(pack_msgpack(payload))
