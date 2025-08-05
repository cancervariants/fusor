#!/usr/bin/env python

import asyncio
from pathlib import Path

from civicpy import civic

from fusor import FUSOR
from fusor.config import config
from fusor.harvester import CIVICHarvester, MOAHarvester
from fusor.models import save_fusions_cache

SYSTEM_CIVIC_CACHE_PATH = Path("~/.civicpy/cache.pkl").expanduser()


async def setup_test_caches():
    ### Set up CIVIC cache
    # Generate list of matches, report match score

    fusor = FUSOR()
    # Generate categorical fusions list
    local_cache_path = (
        SYSTEM_CIVIC_CACHE_PATH if SYSTEM_CIVIC_CACHE_PATH.exists else "civic_cache.pkl"
    )
    harvester = CIVICHarvester(fusor=fusor, local_cache_path=local_cache_path)
    variants = civic.get_all_fusion_variants(include_status="accepted")
    harvester.fusions_list = variants
    civic_fusions = await harvester.load_records()

    # Save cache for later
    save_fusions_cache(
        civic_fusions,
        cache_dir=config.data_root,
        cache_name="civic_translated_fusions.pkl",
    )

    ### Set up MOA cache

    harvester = MOAHarvester(fusor=fusor, cache_dir=Path.cwd())
    moa_fusions = harvester.load_records()

    # Save cache for later
    save_fusions_cache(
        moa_fusions, cache_dir=config.data_root, cache_name="moa_translated_fusions.pkl"
    )


if __name__ == "__main__":
    asyncio.run(setup_test_caches())
