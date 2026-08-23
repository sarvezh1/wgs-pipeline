#!/usr/bin/env python3
"""Capture a local report through a running Chromium DevTools endpoint."""

from __future__ import annotations

import argparse
import base64
import json
import os
import socket
import struct
import urllib.parse
import urllib.request
from pathlib import Path


class WebSocket:
    def __init__(self, url: str):
        parsed = urllib.parse.urlparse(url)
        self.sock = socket.create_connection((parsed.hostname, parsed.port or 80), timeout=15)
        key = base64.b64encode(os.urandom(16)).decode()
        resource = (parsed.path or "/") + (f"?{parsed.query}" if parsed.query else "")
        request = (
            f"GET {resource} HTTP/1.1\r\nHost: {parsed.netloc}\r\n"
            f"Upgrade: websocket\r\nConnection: Upgrade\r\nSec-WebSocket-Key: {key}\r\n"
            "Sec-WebSocket-Version: 13\r\nOrigin: http://localhost\r\n\r\n"
        )
        self.sock.sendall(request.encode())
        response = b""
        while b"\r\n\r\n" not in response:
            response += self.sock.recv(4096)
        if b" 101 " not in response.split(b"\r\n", 1)[0]:
            raise RuntimeError(f"WebSocket upgrade failed: {response[:200]!r}")

    def send(self, message: dict) -> None:
        payload = json.dumps(message, separators=(",", ":")).encode()
        mask = os.urandom(4)
        length = len(payload)
        header = bytearray([0x81])
        if length < 126:
            header.append(0x80 | length)
        elif length < 65536:
            header.extend([0x80 | 126]); header.extend(struct.pack(">H", length))
        else:
            header.extend([0x80 | 127]); header.extend(struct.pack(">Q", length))
        header.extend(mask)
        header.extend(byte ^ mask[index % 4] for index, byte in enumerate(payload))
        self.sock.sendall(header)

    def receive(self) -> dict:
        first, second = self._read(2)
        opcode, length = first & 0x0F, second & 0x7F
        if length == 126:
            length = struct.unpack(">H", self._read(2))[0]
        elif length == 127:
            length = struct.unpack(">Q", self._read(8))[0]
        if second & 0x80:
            mask = self._read(4)
            payload = bytes(byte ^ mask[index % 4] for index, byte in enumerate(self._read(length)))
        else:
            payload = self._read(length)
        if opcode == 8:
            raise RuntimeError("Browser closed the DevTools connection")
        if opcode == 9:
            self._send_control(10, payload)
            return self.receive()
        return json.loads(payload)

    def _read(self, size: int) -> bytes:
        chunks = bytearray()
        while len(chunks) < size:
            chunk = self.sock.recv(size - len(chunks))
            if not chunk:
                raise RuntimeError("Unexpected end of DevTools connection")
            chunks.extend(chunk)
        return bytes(chunks)

    def _send_control(self, opcode: int, payload: bytes) -> None:
        mask = os.urandom(4)
        self.sock.sendall(bytes([0x80 | opcode, 0x80 | len(payload)]) + mask + bytes(byte ^ mask[i % 4] for i, byte in enumerate(payload)))

    def close(self) -> None:
        self.sock.close()


def command(ws: WebSocket, command_id: int, method: str, params: dict | None = None) -> dict:
    ws.send({"id": command_id, "method": method, "params": params or {}})
    while True:
        message = ws.receive()
        if message.get("id") == command_id:
            if "error" in message:
                raise RuntimeError(f"Chromium {method} failed: {message['error']}")
            return message.get("result", {})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint", default="http://127.0.0.1:9222")
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    args = parser.parse_args()
    # Reuse the initial page exposed by headless-shell. Some builds reject
    # /json/new even though their existing DevTools target is fully usable.
    with urllib.request.urlopen(f"{args.endpoint}/json", timeout=15) as response:
        targets = json.load(response)
    target = next((item for item in targets if item.get("type") == "page"), None)
    if target is None:
        request = urllib.request.Request(f"{args.endpoint}/json/new?{urllib.parse.quote(args.url, safe=':/')}", method="PUT")
        with urllib.request.urlopen(request, timeout=15) as response:
            target = json.load(response)
    ws_url = urllib.parse.urlparse(target["webSocketDebuggerUrl"])
    endpoint = urllib.parse.urlparse(args.endpoint)
    ws_url = urllib.parse.urlunparse(ws_url._replace(netloc=endpoint.netloc))
    ws = WebSocket(ws_url)
    try:
        command(ws, 1, "Page.enable")
        command(ws, 2, "Network.enable")
        command(ws, 3, "Emulation.setDeviceMetricsOverride", {"width": args.width, "height": args.height, "deviceScaleFactor": 1, "mobile": False})
        command(ws, 4, "Page.navigate", {"url": args.url})
        while True:
            message = ws.receive()
            if message.get("method") == "Page.loadEventFired":
                break
        result = command(ws, 5, "Page.captureScreenshot", {"format": "png", "fromSurface": True, "captureBeyondViewport": False})
        payload = base64.b64decode(result["data"])
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(payload)
    finally:
        ws.close()
    print(f"Browser-rendered screenshot: {args.output} ({args.width}x{args.height})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
