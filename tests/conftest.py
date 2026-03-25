"""
Pytest configuration and common fixtures
"""

import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing"""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def sample_data():
    """Sample data for testing"""
    return {
        "samples": ["sample1", "sample2", "sample3"],
        "include_list": ["sample1", "sample2"],
        "exclude_list": ["sample3"],
        "file_paths": [
            "/path/to/file1.txt",
            "/path/to/file2.txt",
            "/path/to/dir/file3.txt",
        ],
    }


@pytest.fixture(autouse=True)
def setup_test_environment(monkeypatch):
    """Setup test environment variables and paths"""
    # Set up test-specific environment
    monkeypatch.setenv("TESTING", "true")
