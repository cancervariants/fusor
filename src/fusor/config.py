"""Configure data storage variables for fusion objects"""

import os
from pathlib import Path

from wags_tails.utils.storage import get_data_dir


def get_data_config() -> Path:
    """Get location for storing cached fusion objects

    :return: A Path object
    """
    return (
        Path(env_var_data_dir)
        if (env_var_data_dir := os.environ.get("FUSOR_DIR"))
        else get_data_dir() / "fusor"
    )
