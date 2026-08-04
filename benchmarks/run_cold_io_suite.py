#!/usr/bin/env python3
"""Orchestrate fresh-process, fresh-copy CompreSSoR/FastMR I/O trials."""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
import random
import re
import shutil
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKER = ROOT / "benchmarks" / "cold_io_trial.R"
DEFAULT_INSTRUMENTS = ROOT / "benchmarks" / "finngen_relaxed_instruments.tsv"
MR_FORMATS = ("pcodec_direct", "pcodec_explicit", "vcf_tabix", "tsv_gz")
FULL_FORMATS = ("pcodec_full", "vcf_full", "tsv_full")
SHAPES = {
    "mr_1x1": (1, 1),
    "mr_5x5": (5, 5),
    "mr_25x25": (25, 25),
    "full_load": (0, 0),
}


def open_nocache(path: Path, flags: int, mode: int = 0o600) -> int:
    fd = os.open(path, flags, mode)
    try:
        fcntl.fcntl(fd, fcntl.F_NOCACHE, 1)
    except Exception:
        os.close(fd)
        raise
    return fd


def copy_file_nocache(source: Path, target: Path) -> dict:
    target.parent.mkdir(parents=True, exist_ok=True)
    source_fd = open_nocache(source, os.O_RDONLY)
    target_fd = open_nocache(
        target, os.O_WRONLY | os.O_CREAT | os.O_EXCL, source.stat().st_mode & 0o777
    )
    digest = hashlib.sha256()
    written = 0
    try:
        while True:
            block = os.read(source_fd, 8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
            view = memoryview(block)
            while view:
                count = os.write(target_fd, view)
                if count <= 0:
                    raise OSError(f"short write while staging {target}")
                written += count
                view = view[count:]
        os.fsync(target_fd)
    finally:
        os.close(source_fd)
        os.close(target_fd)
    if written != source.stat().st_size or target.stat().st_size != written:
        raise OSError(f"staged file length differs for {source}")
    return {"bytes": written, "sha256": digest.hexdigest()}


def copy_store_nocache(source: Path, target: Path) -> dict:
    target.mkdir(parents=True, exist_ok=False)
    records = {}
    for path in sorted(source.rglob("*")):
        relative = path.relative_to(source)
        if path.is_dir():
            (target / relative).mkdir(parents=True, exist_ok=True)
        elif path.is_file():
            records[str(relative)] = copy_file_nocache(path, target / relative)
        else:
            raise ValueError(f"unsupported store entry: {path}")
    return records


def regular_bytes(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def load_keys(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = [row["variant_key"] for row in rows]
    if len(keys) != 25 or len(set(keys)) != 25:
        raise ValueError("instrument panel must contain exactly 25 unique keys")
    return keys


def parse_time_output(stderr: str) -> dict:
    fields = {}
    patterns = {
        "maximum_rss_bytes": r"^\s*(\d+)\s+maximum resident set size$",
        "page_faults": r"^\s*(\d+)\s+page faults$",
        "page_reclaims": r"^\s*(\d+)\s+page reclaims$",
        "block_input_operations": r"^\s*(\d+)\s+block input operations$",
        "block_output_operations": r"^\s*(\d+)\s+block output operations$",
    }
    for name, pattern in patterns.items():
        hit = re.search(pattern, stderr, flags=re.MULTILINE)
        fields[name] = int(hit.group(1)) if hit else None
    return fields


def stage_paths(args, input_format: str, count: int,
                trial: Path) -> tuple[list[Path], int]:
    paths = []
    total = 0
    for index in range(count):
        label = f"study-{index + 1:02d}"
        if input_format.startswith("pcodec"):
            target = args.store
        elif input_format.startswith("tsv"):
            target = trial / f"{label}.tsv.gz"
            record = copy_file_nocache(args.tsv, target)
            total += int(record["bytes"])
        elif input_format.startswith("vcf"):
            target = trial / f"{label}.vcf.gz"
            record = copy_file_nocache(args.vcf, target)
            index_record = copy_file_nocache(
                Path(str(args.vcf) + ".tbi"), Path(str(target) + ".tbi")
            )
            total += int(record["bytes"]) + int(index_record["bytes"])
        else:
            raise ValueError(f"unknown format: {input_format}")
        paths.append(target)
    return paths, total


def run_trial(args, keys: list[str], repetition: int, workload: str,
              input_format: str, order: int) -> dict:
    exposures, outcomes = SHAPES[workload]
    study_count = 1 if workload == "full_load" else exposures + outcomes
    trial = Path(tempfile.mkdtemp(
        prefix=f"{repetition:02d}-{workload}-{input_format}-", dir=args.scratch
    ))
    bridge_temp = Path(tempfile.mkdtemp(prefix="compressor-cold-bridge-"))
    started = time.perf_counter()
    try:
        paths, staged_bytes = stage_paths(args, input_format, study_count, trial)
        staging_seconds = time.perf_counter() - started
        spec = {
            "format": input_format,
            "workload": workload,
            "paths": [str(path) for path in paths],
            "keys": keys,
            "exposures": exposures,
            "outcomes": outcomes,
            "io_threads": (
                min(args.pcodec_io_threads, study_count)
                if input_format == "pcodec_direct" else 1
            ),
            "tabix": str(args.tabix),
            "bcftools": str(args.bcftools),
            "tempdir": str(bridge_temp),
        }
        spec_path = trial / "spec.json"
        result_path = trial / "result.json"
        spec_path.write_text(json.dumps(spec, indent=2) + "\n")
        environment = os.environ.copy()
        environment.update({
            "R_LIBS_USER": str(args.r_library),
            "COMPRESSOR_PYTHON": str(args.python),
            "COMPRESSOR_NOCACHE": "1",
            "COMPRESSOR_PAGE_CACHE_PAGES": "0",
            "COMPRESSOR_TMPDIR": str(bridge_temp),
            "TMPDIR": str(bridge_temp),
        })
        command = [
            "/usr/bin/time", "-l", str(args.rscript), str(WORKER),
            str(spec_path), str(result_path),
        ]
        process = subprocess.run(
            command, env=environment, text=True, capture_output=True
        )
        if process.returncode:
            raise RuntimeError(
                f"trial failed ({workload}, {input_format}):\n"
                f"STDOUT:\n{process.stdout}\nSTDERR:\n{process.stderr}"
            )
        record = json.loads(result_path.read_text())
        record.update(parse_time_output(process.stderr))
        record.update({
            "repetition": repetition,
            "randomized_order": order,
            "exposures": exposures,
            "outcomes": outcomes,
            "logical_study_reads": study_count,
            "io_threads": int(record.get("io_threads", 1)),
            # Staging is deliberately outside the timed interval. Pcodec uses
            # F_NOCACHE on the immutable source itself, so its staging is zero.
            "physical_staging_bytes_not_timed": staged_bytes,
            "staging_seconds_not_timed": staging_seconds,
            "cache_protocol": (
                "fresh R process; every logical Pcodec study gets a fresh reader "
                "context with coalescing/page cache disabled, exceptions reloaded, "
                "and F_NOCACHE; TSV/VCF use fresh F_NOCACHE-staged copies read once"
            ),
        })
        validate_record(record, keys)
        return record
    finally:
        shutil.rmtree(trial, ignore_errors=True)
        shutil.rmtree(bridge_temp, ignore_errors=True)


def mr_rows(record: dict) -> list[dict]:
    value = record["mr"]
    if isinstance(value, list):
        return value
    if isinstance(value, dict):
        names = list(value)
        length = len(value[names[0]])
        return [{name: value[name][index] for name in names} for index in range(length)]
    raise ValueError("MR result has an unsupported JSON shape")


def validate_record(record: dict, keys: list[str]) -> None:
    if record["workload"] == "full_load":
        checksum = record["checksum"]
        if int(checksum["rows"]) != 14_923_434:
            raise ValueError("full-load row count differs from the release fixture")
        if not checksum["first_key"] or not checksum["last_key"]:
            raise ValueError("full-load identity bounds are empty")
        return
    rows = mr_rows(record)
    expected_pairs = int(record["exposures"]) * int(record["outcomes"])
    if len(rows) != expected_pairs:
        raise ValueError("MR result has the wrong number of exposure/outcome pairs")
    for row in rows:
        if int(row["nsnp"]) != len(keys):
            raise ValueError("MR pair did not retain all instruments")
        if not all(math.isfinite(float(row[name])) for name in (
                "b", "se", "Q", "Q_pval", "sigma")):
            raise ValueError("MR result contains non-finite estimate or standard error")
        if not math.isclose(float(row["b"]), 1.0, rel_tol=0, abs_tol=1e-12):
            raise ValueError("self-GWAS MR estimate differs from one")
        if not math.isclose(float(row["se"]), 0.0, rel_tol=0, abs_tol=1e-12):
            raise ValueError(
                f"self-GWAS MR standard error differs from zero: {row!r}"
            )


def mr_signature(record: dict) -> list[tuple]:
    return [
        (
            row["id_exposure"], row["id_outcome"], row["method"],
            int(row["nsnp"]), float(row["b"]), float(row["se"]),
            # jsonlite omits non-finite data-frame fields.  A self-GWAS has
            # b=1 and se=0, so its p-value is intentionally undefined.
            None if row.get("pval") is None else float(row["pval"]),
            float(row["Q"]), float(row["Q_pval"]), float(row["sigma"]),
        )
        for row in mr_rows(record)
    ]


def parity_report(records: list[dict]) -> list[dict]:
    report = []
    for workload in sorted({row["workload"] for row in records}):
        selected = [row for row in records if row["workload"] == workload]
        if workload == "full_load":
            exact = [row for row in selected if row["format"] in ("tsv_full", "vcf_full")]
            identity_values = {
                (
                    int(row["checksum"]["rows"]),
                    float(row["checksum"]["identity_sum"]),
                    row["checksum"]["first_key"], row["checksum"]["last_key"],
                )
                for row in selected
            }
            exact_numeric = {
                tuple(float(row["checksum"][name]) for name in (
                    "beta_sum", "beta_sum_squares", "se_sum", "eaf_sum", "z_sum"
                ))
                for row in exact
            }
            missing_counts = {
                tuple(int(row["checksum"][name]) for name in (
                    "beta_missing", "se_missing", "eaf_missing", "z_missing"
                ))
                for row in selected
            }
            report.append({
                "workload": workload,
                "identity_exact_all_formats": len(identity_values) == 1,
                "tsv_vcf_numeric_exact": len(exact_numeric) <= 1,
                "missing_counts_exact_all_formats": len(missing_counts) == 1,
            })
            continue
        groups = {}
        format_records = {}
        for row in selected:
            groups.setdefault(row["format"], []).append(mr_signature(row))
            format_records.setdefault(row["format"], []).append(row)
        pcodec_exact = len({repr(value) for value in (
            groups.get("pcodec_direct", []) + groups.get("pcodec_explicit", [])
        )}) <= 1
        exact_formats = groups.get("tsv_gz", []) + groups.get("vcf_tabix", [])
        exact_equal = len({repr(value) for value in exact_formats}) <= 1
        signatures = {
            name: values[0]["input_signature"]
            for name, values in format_records.items()
            if values and "input_signature" in values[0]
        }
        raw_signature = signatures.get("tsv_gz")
        vcf_signature = signatures.get("vcf_tabix")
        pcodec_signature = signatures.get("pcodec_explicit")
        raw_exact = (
            raw_signature is not None and vcf_signature is not None and
            raw_signature == vcf_signature
        )
        profile_ok = False
        numeric_delta = None
        if raw_signature is not None and pcodec_signature is not None:
            raw_keys = list(raw_signature["variant_key"])
            pcodec_keys = list(pcodec_signature["variant_key"])
            if raw_keys == pcodec_keys:
                raw_beta = [float(value) for value in raw_signature["beta"]]
                raw_se = [float(value) for value in raw_signature["standard_error"]]
                packed_beta = [float(value) for value in pcodec_signature["beta"]]
                packed_se = [
                    float(value) for value in pcodec_signature["standard_error"]
                ]
                beta_delta = max(abs(a - b) for a, b in zip(raw_beta, packed_beta))
                se_delta = max(abs(a - b) for a, b in zip(raw_se, packed_se))
                z_delta = max(abs(
                    a / sa - b / sb
                ) for a, sa, b, sb in zip(
                    raw_beta, raw_se, packed_beta, packed_se
                ))
                numeric_delta = {
                    "max_abs_beta": beta_delta,
                    "max_abs_standard_error": se_delta,
                    "max_abs_z": z_delta,
                }
                # Z9 has a central half-bin error of 7/(2*510); exception
                # values are float32. The broad beta/SE gates are the package's
                # documented standard-profile regression tolerances.
                profile_ok = (
                    z_delta <= 7.0 / (2.0 * 510.0) + 1e-6 and
                    beta_delta <= 0.02 and se_delta <= 0.01
                )
        entry = {
            "workload": workload,
            "pcodec_direct_explicit_exact": pcodec_exact,
            "tsv_vcf_exact": exact_equal,
        }
        if raw_signature is not None and vcf_signature is not None:
            entry["tsv_vcf_input_exact"] = raw_exact
        if raw_signature is not None and pcodec_signature is not None:
            entry["pcodec_raw_numeric_within_profile"] = profile_ok
            entry["pcodec_raw_max_delta"] = numeric_delta
        report.append(entry)
    return report


def summaries(records: list[dict]) -> list[dict]:
    def optional_number(value):
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            return None
        return parsed if math.isfinite(parsed) else None

    groups = {}
    for record in records:
        groups.setdefault((record["workload"], record["format"]), []).append(record)
    output = []
    for (workload, input_format), values in sorted(groups.items()):
        total = [float(row["total_seconds"]) for row in values]
        io = [float(row["io_seconds"]) for row in values]
        estimator = [float(row["estimator_seconds"]) for row in values]
        rss = [int(row["maximum_rss_bytes"]) for row in values]
        source_read = [
            parsed for parsed in (
                optional_number(row.get("source_bytes_read")) for row in values
            ) if parsed is not None
        ]
        output.append({
            "workload": workload,
            "format": input_format,
            "runs": len(values),
            "io_threads": int(values[0].get("io_threads", 1)),
            "median_total_seconds": statistics.median(total),
            "min_total_seconds": min(total),
            "max_total_seconds": max(total),
            "median_io_seconds": statistics.median(io),
            "median_estimator_seconds": statistics.median(estimator),
            "median_peak_rss_bytes": statistics.median(rss),
            "median_pcodec_source_bytes_read": (
                statistics.median(source_read) if source_read else None
            ),
        })
    return output


def git_revision(path: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument("--tsv", type=Path, required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--tabix", type=Path, default=Path("/opt/homebrew/bin/tabix"))
    parser.add_argument("--bcftools", type=Path, default=Path("/opt/homebrew/bin/bcftools"))
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--rscript", type=Path, default=Path("/opt/homebrew/bin/Rscript"))
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument(
        "--compressor-repo", type=Path, default=ROOT.parent / "CompreSSoR"
    )
    parser.add_argument("--instruments", type=Path, default=DEFAULT_INSTRUMENTS)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument(
        "--pcodec-io-threads", type=int, default=1,
        help="parallel store decoders used by the direct Pcodec/FastMR path"
    )
    parser.add_argument("--workload", action="append", choices=tuple(SHAPES))
    parser.add_argument(
        "--format", action="append", dest="formats",
        choices=MR_FORMATS + FULL_FORMATS,
        help="restrict diagnostic runs to selected compatible format(s)"
    )
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.reps < 5 and not args.smoke:
        raise ValueError("release benchmark requires at least five repetitions")
    if args.pcodec_io_threads < 1:
        raise ValueError("--pcodec-io-threads must be positive")
    required = [
        args.store, args.tsv, args.vcf, Path(str(args.vcf) + ".tbi"),
        args.tabix, args.bcftools, args.python, args.rscript, args.r_library,
        args.instruments, args.compressor_repo, WORKER,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing benchmark input(s): " + ", ".join(missing))
    args.scratch.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    keys = load_keys(args.instruments)
    workloads = args.workload or list(SHAPES)
    records = []
    schedule = []
    rng = random.Random(args.seed)
    for repetition in range(1, args.reps + 1):
        conditions = []
        for workload in workloads:
            formats = list(FULL_FORMATS if workload == "full_load" else MR_FORMATS)
            if args.formats:
                formats = [value for value in formats if value in args.formats]
                if not formats:
                    raise ValueError(
                        f"none of the selected formats is compatible with {workload}"
                    )
            conditions.extend((workload, value) for value in formats)
        rng.shuffle(conditions)
        schedule.extend(
            (repetition, workload, value) for workload, value in conditions
        )
    for order, (repetition, workload, input_format) in enumerate(schedule, 1):
        print(
            f"trial {order}/{len(schedule)}: rep={repetition} "
            f"workload={workload} format={input_format}", flush=True
        )
        record = run_trial(
            args, keys, repetition, workload, input_format, order
        )
        records.append(record)
        print(
            f"completed {workload} {input_format}: "
            f"{record['total_seconds']:.6f}s", flush=True
        )
        checkpoint = {
            "schema_version": "1.0.0",
            "complete": False,
            "records": records,
        }
        (args.output / "cold-io-runs.partial.json").write_text(
            json.dumps(checkpoint, indent=2) + "\n"
        )

    parity = parity_report(records)
    gate_suffixes = ("_exact", "_exact_all_formats", "_within_profile")
    if not all(all(
            value is True for key, value in row.items()
            if key.endswith(gate_suffixes)
    ) for row in parity):
        raise RuntimeError(f"parity gate failed: {parity}")
    summary = summaries(records)
    result = {
        "schema_version": "1.0.0",
        "complete": True,
        "measured_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "seed": args.seed,
        "repetitions": args.reps,
        "pcodec_io_threads": args.pcodec_io_threads,
        "instrument_panel": {
            "path": str(args.instruments),
            "sha256": hashlib.sha256(args.instruments.read_bytes()).hexdigest(),
            "variants": len(keys),
            "selection": (
                "FinnGen p<1e-5; PLINK2 clump against GRCh38 1KG EUR, "
                "r2<0.001, 10000 kb; first 25 index variants by p"
            ),
        },
        "source_revisions": {
            "CompreSSoR": git_revision(args.compressor_repo),
            "fastMR": git_revision(ROOT),
        },
        "source_bytes": {
            "pcodec": regular_bytes(args.store),
            "tsv_gz": regular_bytes(args.tsv),
            "vcf_gz_plus_tbi": regular_bytes(args.vcf) + regular_bytes(
                Path(str(args.vcf) + ".tbi")
            ),
        },
        "cache_protocol": (
            "fresh R process; every logical Pcodec study gets a fresh reader "
            "context for the same immutable store, with request coalescing/page "
            "cache disabled, exceptions reloaded, and F_NOCACHE on data "
            "descriptors; TSV/VCF use fresh copies written through F_NOCACHE "
            "descriptors and read once (a cache-controlled cold approximation)"
        ),
        "parity": parity,
        "summary": summary,
        "records": records,
    }
    (args.output / "cold-io-benchmark.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    with (args.output / "cold-io-summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    with (args.output / "cold-io-runs.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = [
            "repetition", "randomized_order", "workload", "format",
            "exposures", "outcomes", "logical_study_reads", "io_threads",
            "total_seconds",
            "io_seconds", "estimator_seconds", "touch_seconds",
            "source_bytes_read",
            "maximum_rss_bytes", "page_faults", "page_reclaims",
            "block_input_operations", "block_output_operations",
            "physical_staging_bytes_not_timed",
            "staging_seconds_not_timed",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    partial = args.output / "cold-io-runs.partial.json"
    if partial.exists():
        partial.unlink()
    print(json.dumps({"summary": summary, "parity": parity}, indent=2))


if __name__ == "__main__":
    main()
