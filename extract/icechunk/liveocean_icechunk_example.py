"""
Create an example virtual Icechunk store for LiveOcean forecast files.

This intentionally avoids listing the source bucket. It derives daily
``layers.nc`` paths from a known run path and processes the seven preceding
runs by default.

Originally from Rich Signell, 5/2026. Modified by PM to use kopeh buckets for output.

modules added to loenv:
icechunk, obspec_utils, obstore, virtualizarr
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import datetime as dt
import os
import re
import warnings
from pathlib import PurePosixPath

import icechunk
import xarray as xr
from dotenv import load_dotenv
from obspec_utils.registry import ObjectStoreRegistry
from obstore.store import S3Store
from virtualizarr import open_virtual_dataset
from virtualizarr.parsers import HDFParser

warnings.filterwarnings("ignore", category=UserWarning)

KNOWN_SOURCE_URL = "s3://liveocean-share/f2026.05.08/layers.nc"
SOURCE_BUCKET = "liveocean-share"
SOURCE_BUCKET_URL = f"s3://{SOURCE_BUCKET}"
SOURCE_ENDPOINT_URL = "https://s3.kopah.uw.edu"
TARGET_STORAGE_NAME = "liveocean-layers-icechunk-example"


@dataclass(frozen=True)
class TargetConfig:
    bucket: str
    prefix: str
    endpoint_url: str | None
    region: str
    force_path_style: bool


def parse_run_date(url: str) -> dt.date:
    match = re.search(r"/f(\d{4})\.(\d{2})\.(\d{2})/", url)
    if not match:
        raise ValueError(f"Could not find /fYYYY.MM.DD/ run date in {url!r}")

    year, month, day = (int(part) for part in match.groups())
    return dt.date(year, month, day)


def build_source_urls(known_url: str, days_before: int = 7) -> list[str]:
    run_date = parse_run_date(known_url)
    source_path = PurePosixPath(known_url.removeprefix(f"{SOURCE_BUCKET_URL}/"))
    filename = source_path.name

    urls = []
    for offset in range(days_before, 0, -1):
        date = run_date - dt.timedelta(days=offset)
        urls.append(f"{SOURCE_BUCKET_URL}/f{date:%Y.%m.%d}/{filename}")
    return urls


def standardize_forecast_run(ds: xr.Dataset) -> xr.Dataset:
    """Convert LiveOcean ocean_time to forecast reference time plus step."""
    valid_time = ds["ocean_time"]
    step = valid_time - valid_time[0]
    time = valid_time[0]

    ds = ds.rename_dims({"ocean_time": "step"})
    ds = ds.drop_vars("ocean_time")
    ds = ds.assign_coords(
        step=("step", step.data, {"standard_name": "forecast_period"}),
        time=time.assign_attrs({"standard_name": "forecast_reference_time"}),
    )
    return ds


def promote_lon_lat_coordinates(ds: xr.Dataset) -> xr.Dataset:
    coord_names = [
        name
        for name, data_array in ds.data_vars.items()
        if (name.startswith("lon_") or name.startswith("lat_"))
        and "time" not in data_array.dims
        and "step" not in data_array.dims
    ]
    if not coord_names:
        return ds
    return ds.assign_coords({name: ds[name] for name in coord_names})


def concat_forecast_runs(datasets: list[xr.Dataset]) -> xr.Dataset:
    datasets = [promote_lon_lat_coordinates(ds) for ds in datasets]
    return xr.concat(
        datasets,
        dim="time",
        coords="minimal",
        compat="override",
        combine_attrs="override",
    )


def make_source_registry() -> ObjectStoreRegistry:
    source_store = S3Store(
        bucket=SOURCE_BUCKET,
        endpoint=SOURCE_ENDPOINT_URL,
        region="not-used",
        skip_signature=True,
    )
    return ObjectStoreRegistry({SOURCE_BUCKET_URL: source_store})


def make_target_s3_credentials():
    access_key_id = os.environ.get("AWS_ACCESS_KEY_ID")
    secret_access_key = os.environ.get("AWS_SECRET_ACCESS_KEY")
    # session_token = os.environ.get("AWS_SESSION_TOKEN")

    missing = [
        name
        for name, value in (
            ("AWS_ACCESS_KEY_ID", access_key_id),
            ("AWS_SECRET_ACCESS_KEY", secret_access_key),
        )
        if not value
    ]
    if missing:
        raise RuntimeError(
            "Missing target S3 write credential(s) from "
            f"the selected env file: {', '.join(missing)}"
        )

    return icechunk.s3_credentials(
        access_key_id=access_key_id,
        secret_access_key=secret_access_key,
        # session_token=session_token,
    )


def require_env(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"Missing required environment variable: {name}")
    return value


def get_target_config(target: str, storage_name: str) -> TargetConfig:
    if target == "source-coop":
        return TargetConfig(
            bucket=require_env("SOURCE_COOP_BUCKET"),
            prefix=require_env("SOURCE_COOP_PREFIX").strip("/"),
            endpoint_url=None,
            region=os.environ.get("AWS_DEFAULT_REGION", "us-west-2"),
            force_path_style=False,
        )

    if target == "liveocean-source":
        return TargetConfig(
            bucket=SOURCE_BUCKET,
            prefix=require_env("LIVEOCEAN_ICECHUNK_PREFIX").strip("/"),
            endpoint_url=os.environ.get(
                "LIVEOCEAN_S3_ENDPOINT_URL",
                SOURCE_ENDPOINT_URL,
            ),
            region=os.environ.get("LIVEOCEAN_S3_REGION", "not-used"),
            force_path_style=True,
        )

    raise ValueError(f"Unknown target {target!r}")


def make_icechunk_storage(
    target_config: TargetConfig,
) -> tuple[icechunk.Storage, icechunk.RepositoryConfig, dict]:
    target_bucket_url = f"s3://{target_config.bucket}"
    storage = icechunk.s3_storage(
        bucket=target_config.bucket,
        prefix=target_config.prefix,
        from_env=True,
        endpoint_url=target_config.endpoint_url,
        region=target_config.region,
        force_path_style=target_config.force_path_style,
    )

    config = icechunk.RepositoryConfig.default()
    config.set_virtual_chunk_container(
        icechunk.VirtualChunkContainer(
            url_prefix=f"{SOURCE_BUCKET_URL}/",
            store=icechunk.s3_store(
                region="not-used",
                anonymous=True,
                s3_compatible=True,
                force_path_style=True,
                endpoint_url=SOURCE_ENDPOINT_URL,
            ),
        )
    )

    credentials = icechunk.containers_credentials(
        {
            f"{target_bucket_url}/": make_target_s3_credentials(),
            f"{SOURCE_BUCKET_URL}/": icechunk.s3_credentials(anonymous=True),
        }
    )
    return storage, config, credentials


def build_virtual_dataset(urls: list[str]) -> xr.Dataset:
    parser = HDFParser()
    registry = make_source_registry()
    datasets = [
        standardize_forecast_run(
            open_virtual_dataset(
                url,
                parser=parser,
                registry=registry,
                loadable_variables=["ocean_time"],
            )
        )
        for url in urls
    ]
    return concat_forecast_runs(datasets)


def write_icechunk(ds: xr.Dataset, storage_name: str, mode: str) -> None:
    target_config = get_target_config(os.environ["LIVE_OCEAN_TARGET"], storage_name)
    storage, config, credentials = make_icechunk_storage(target_config)

    if mode == "create":
        repo = icechunk.Repository.create(
            storage,
            config,
            authorize_virtual_chunk_access=credentials,
        )
        append_dim = None
    else:
        repo = icechunk.Repository.open(
            storage,
            config,
            authorize_virtual_chunk_access=credentials,
        )
        append_dim = "time"

    session = repo.writable_session("main")
    if append_dim is None:
        ds.virtualize.to_icechunk(session.store)
    else:
        ds.virtualize.to_icechunk(session.store, append_dim=append_dim)

    first = str(ds.time.values[0])[:10]
    last = str(ds.time.values[-1])[:10]
    session.commit(f"LiveOcean layers example {mode}: {first} to {last}")

    latest = next(repo.ancestry(branch="main"))
    print(f"Committed [{latest.written_at}]: {latest.message}")
    print(f"Store: s3://{target_config.bucket}/{target_config.prefix}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--known-url", default=KNOWN_SOURCE_URL)
    parser.add_argument("--days-before", type=int, default=7)
    parser.add_argument("--storage-name", default=TARGET_STORAGE_NAME)
    parser.add_argument("--mode", choices=["create", "append"], default="create")
    parser.add_argument(
        "--target",
        choices=["source-coop", "liveocean-source"],
        default="source-coop",
        help="Destination object store for the Icechunk metadata.",
    )
    parser.add_argument(
        "--env-file",
        default=None,
        help="Env file with destination credentials. Defaults depend on --target.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the derived URLs without reading source files or writing Icechunk.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    default_env_files = {
        "source-coop": f'{os.environ["HOME"]}/dotenv/source-coop-liveocean.env',
        "liveocean-source": f'{os.environ["HOME"]}/dotenv/liveocean-source.env',
    }
    load_dotenv(args.env_file or default_env_files[args.target], override=True)
    os.environ["LIVE_OCEAN_TARGET"] = args.target

    urls = build_source_urls(args.known_url, days_before=args.days_before)
    print("LiveOcean source URLs:")
    for url in urls:
        print(f"  {url}")

    if args.dry_run:
        return

    ds = build_virtual_dataset(urls)
    print(ds)
    write_icechunk(ds, storage_name=args.storage_name, mode=args.mode)


if __name__ == "__main__":
    main()