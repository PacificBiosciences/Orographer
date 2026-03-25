#!/usr/bin/env python
"""
Deploy module for serving Bokeh plots with external JSON files.

This module provides functionality to start a simple HTTP server that serves
generated Bokeh plots.
"""

import http.server
import logging
import os
import socketserver
import sys

logger = logging.getLogger(__name__)


class PlotHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    """Simple HTTP handler for serving plot files."""

    def end_headers(self):
        # Add CORS headers to allow local file access
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        super().end_headers()


def run_deploy(output_dir, port=8000):
    """
    Start HTTP server to serve generated plots.

    Args:
        output_dir (str): Directory path containing HTML and JSON files to serve
        port (int): Port number to serve on (default: 8000)
    """
    directory = os.path.abspath(output_dir)

    if not os.path.isdir(directory):
        logger.error(f"Error: Directory not found: {directory}")
        sys.exit(1)

    # Change to the directory to serve
    os.chdir(directory)

    # Start server with simple HTTP handler
    with socketserver.TCPServer(("", port), PlotHTTPRequestHandler) as httpd:
        print(f"Serving plots from: {directory}")
        print(f"Server running at http://localhost:{port}/")
        print("\nPress Ctrl+C to stop the server")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n\nServer stopped.")
