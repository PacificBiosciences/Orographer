"""Static lazy-load store assembly for IsoSeq plot data."""

import logging
import os
import time

from . import coverage as coverage_data
from . import isoseq as isoseq_data
from .lazy_store import write_json_gzip, write_msgpack_gzip

logger = logging.getLogger(__name__)


def read_manifest_url(chunk_url_prefix: str | None, chunk_key: str) -> str:
    """Return a browser-loadable static manifest URL for IsoSeq read pages."""
    manifest_file = f"{chunk_key}_manifest.json.gz"
    return f"{chunk_url_prefix}/{manifest_file}" if chunk_url_prefix else ""


def all_assignments_url(chunk_url_prefix: str | None, assignments_key: str) -> str:
    """Return a browser-loadable URL for the cross-annotation read assignments file."""
    return f"{chunk_url_prefix}/{assignments_key}_all_assignments.json.gz" if chunk_url_prefix else ""


def write_all_assignments(
    chunk_dir: str | None,
    assignments_key: str,
    all_read_assignments: dict[str, list[dict]],
) -> None:
    """Write the cross-annotation read-to-isoform assignment lookup file."""
    if not chunk_dir or not all_read_assignments:
        return
    os.makedirs(chunk_dir, exist_ok=True)
    path = os.path.join(chunk_dir, f"{assignments_key}_all_assignments.json.gz")
    write_json_gzip(path, {"assignments": all_read_assignments})


def coverage_url(chunk_url_prefix: str | None, chunk_key: str) -> str:
    """Return a browser-loadable static JSON URL for one transcript coverage payload."""
    return f"{chunk_url_prefix}/{chunk_key}_coverage.json.gz" if chunk_url_prefix else ""


def write_coverage_page(
    chunk_dir: str | None,
    chunk_key: str,
    coverage_source_data: dict,
) -> None:
    """Write one static IsoSeq transcript coverage payload."""
    if not chunk_dir:
        return
    os.makedirs(chunk_dir, exist_ok=True)
    coverage_path = os.path.join(chunk_dir, f"{chunk_key}_coverage.json.gz")
    write_json_gzip(coverage_path, {"coverage_data": coverage_source_data})


def write_read_store(
    chunk_dir: str | None,
    manifest_key: str,
    groups: list[dict],
    read_records: list,
    read_id_by_name: dict[str, int],
    manifest_context: dict,
    read_page_size: int = isoseq_data.ISOSEQ_READ_PAGE_SIZE,
) -> None:
    """Write sharded read records plus the read manifest for one BAM display row."""
    if not chunk_dir:
        return
    os.makedirs(chunk_dir, exist_ok=True)
    shard_files = []
    shard_size = isoseq_data.ISOSEQ_READ_SHARD_SIZE
    for shard_index in range(0, len(read_records), shard_size):
        shard_id = shard_index // shard_size
        shard_file = f"{manifest_key}_shard{shard_id}.msgpack.gz"
        shard_path = os.path.join(chunk_dir, shard_file)
        shard_records = read_records[shard_index : shard_index + shard_size]
        write_msgpack_gzip(
            shard_path,
            {
                "schema": "isoseq_read_shard_v1",
                "shard_id": shard_id,
                "reads": shard_records,
            },
        )
        shard_files.append(shard_file)

    manifest_groups = {}
    for group in groups:
        chunk_key = group.get("chunk_key", "")
        if not chunk_key:
            continue
        manifest_groups[chunk_key] = {
            "read_ids": [
                read_id_by_name[name]
                for name in group.get("read_names", [])
                if name in read_id_by_name
            ],
            "assigned_read_count": int(group.get("assigned_read_count", 0)),
        }
    manifest = {
        "schema": "isoseq_read_manifest_v1",
        "page_size": int(read_page_size),
        "shard_size": shard_size,
        "shards": shard_files,
        "groups": manifest_groups,
        **manifest_context,
    }
    manifest_path = os.path.join(chunk_dir, f"{manifest_key}_manifest.json.gz")
    write_json_gzip(manifest_path, manifest)


def prepare_lazy_chunks(
    groups: list[dict],
    segments_by_read: dict,
    layout: dict,
    coordinate_start: int,
    coordinate_end: int,
    region_idx: int,
    row_index: int,
    sample_label: str | None,
    annotation_id: str = "primary",
    annotation_label: str = "Primary GTF",
    assignments_key: str | None = None,
    chunk_dir: str | None = None,
    chunk_url_prefix: str | None = None,
    read_page_size: int = isoseq_data.ISOSEQ_READ_PAGE_SIZE,
) -> dict:
    """Write static sharded read data for browser-side IsoSeq lazy loading."""
    chunk_start_time = time.perf_counter()
    coverage_blocks = coverage_data.coverage_block_cache(
        segments_by_read,
        coordinate_start,
        coordinate_end,
    )
    total_coverage_data = coverage_data.coverage_for_cached_reads(
        list(segments_by_read),
        coverage_blocks,
        coordinate_start,
        coordinate_end,
    )
    if annotation_id == "primary":
        manifest_key = f"r{region_idx}_row{row_index}_reads"
    else:
        annotation_token = isoseq_data.safe_chunk_token(annotation_id)
        manifest_key = f"r{region_idx}_row{row_index}_{annotation_token}_reads"
    manifest_url = read_manifest_url(chunk_url_prefix, manifest_key)
    read_records = []
    read_id_by_name = {}
    for group_index, group in enumerate(groups):
        assigned_read_names = [
            name for name in group.get("read_names", []) if name in segments_by_read
        ]
        chunk_key = f"g{group_index}_{isoseq_data.safe_chunk_token(group.get('group_id'))}"
        group["chunk_key"] = chunk_key
        if not assigned_read_names:
            group["chunk_url"] = ""
            group["coverage_url"] = ""
            continue
        read_metadata = {
            read_name: layout["read_metadata"].get(read_name, {})
            for read_name in assigned_read_names
        }
        for read_name in assigned_read_names:
            if read_name in read_id_by_name:
                continue
            read_id_by_name[read_name] = len(read_records)
            read_records.append(
                isoseq_data.read_shard_record(
                    read_name,
                    segments_by_read[read_name],
                    read_metadata.get(read_name, {}),
                )
            )
        coverage_source_data = coverage_data.coverage_for_cached_reads(
            assigned_read_names,
            coverage_blocks,
            coordinate_start,
            coordinate_end,
        )
        write_coverage_page(chunk_dir, chunk_key, coverage_source_data)
        group["chunk_url"] = manifest_url
        group["coverage_url"] = coverage_url(chunk_url_prefix, chunk_key)
    if read_records:
        write_read_store(
            chunk_dir,
            manifest_key,
            groups,
            read_records,
            read_id_by_name,
            {
                "chrom": str(next(iter(segments_by_read.values()))[0].chrom),
                "coordinate_start": int(coordinate_start),
                "coordinate_end": int(coordinate_end),
                "region_idx": int(region_idx),
                "row_index": int(row_index),
                "sample_label": sample_label or "",
                "annotation_label": annotation_label,
                "all_assignments_url": all_assignments_url(chunk_url_prefix, assignments_key)
                if assignments_key else "",
                "selected_read_y_start": layout["selected_read_y_start"],
            },
            read_page_size=read_page_size,
        )
    logger.debug(
        "IsoSeq static read-store generation for region %d row %d wrote %d read(s) in %.3fs.",
        region_idx,
        row_index,
        len(read_records),
        time.perf_counter() - chunk_start_time,
    )
    return total_coverage_data
