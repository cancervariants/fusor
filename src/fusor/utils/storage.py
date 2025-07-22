"""Provide function for configuring FUSOR data"""

import os
from pathlib import Path


def get_data_dir() -> Path:
    """Configure location of cached FUSOR data (e.g. categorical fusion objects).

    By default, conform to `XDG Base Directory Specification <https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html>`_,
    unless a directory is specified otherwise:

    1) check env var ``"FUSOR_DIR"``
    2) check env var ``"XDG_DATA_HOME"``
    3) check env var ``"XDG_DATA_DIRS"`` for a colon-separated list, skipping any
        that can't be used (i.e. they're already a file)
    4) otherwise, use ``~/.local/share/``

    :return: path to base data directory
    """
    spec_fusor_dir = os.environ.get("FUSOR_DIR")
    if spec_fusor_dir:
        data_base_dir = Path(spec_fusor_dir)
    else:
        xdg_data_home = os.environ.get("XDG_DATA_HOME")
        if xdg_data_home:
            data_base_dir = Path(xdg_data_home) / "fusor"
        else:
            xdg_data_dirs = os.environ.get("XDG_DATA_DIRS")
            if xdg_data_dirs:
                dirs = os.environ["XDG_DATA_DIRS"].split(":")
                for directory in dirs:
                    dir_path = Path(directory) / "fusor"
                    if not dir_path.is_file():
                        data_base_dir = dir_path
                        break
                else:
                    data_base_dir = Path.home() / ".local" / "share" / "fusor"
            else:
                data_base_dir = Path.home() / ".local" / "share" / "fusor"

    data_base_dir.mkdir(exist_ok=True, parents=True)
    return data_base_dir
