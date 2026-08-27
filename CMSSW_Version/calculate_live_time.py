#!/usr/bin/env python3
"""Calculate real-data exposure from muon-tomography ROOT output.

The script sums the RunMetadata equivalent live time from one or more ROOT
files and counts a GEM coincidence directly from MuonHits.  It is intended for
files merged with hadd, where RunMetadata has one row per original job.

Examples:
  python3 calculate_live_time.py MoundTomographyData_AirRoom.root
  python3 calculate_live_time.py *.root --selection five-gem --target-events 100000
  python3 calculate_live_time.py data.root --branches GEMIn2_Valid GEMOut1_Valid
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import uproot


SELECTIONS = {
    "in2-out123": (
        "GEMIn2_Valid, GEMOut1_Valid, GEMOut2_Valid, GEMOut3_Valid",
        ("GEMIn2_Valid", "GEMOut1_Valid", "GEMOut2_Valid", "GEMOut3_Valid"),
    ),
    "five-gem": (
        "GEMIn1_Valid, GEMIn2_Valid, GEMOut1_Valid, GEMOut2_Valid, GEMOut3_Valid",
        ("GEMIn1_Valid", "GEMIn2_Valid", "GEMOut1_Valid", "GEMOut2_Valid", "GEMOut3_Valid"),
    ),
}


def metadata_live_time(metadata: uproot.behaviors.TTree.TTree) -> tuple[int, float]:
    """Return number of merged jobs and summed equivalent live time in seconds."""
    if "EquivalentLiveTime_s" not in metadata:
        raise ValueError("RunMetadata has no EquivalentLiveTime_s branch")
    live_time = metadata["EquivalentLiveTime_s"].array(library="np")
    return int(metadata.num_entries), float(np.sum(live_time))


def count_coincidences(hits: uproot.behaviors.TTree.TTree,
                       branches: tuple[str, ...], step_size: str) -> int:
    missing = [branch for branch in branches if branch not in hits]
    if missing:
        raise ValueError("MuonHits is missing branch(es): " + ", ".join(missing))

    selected = 0
    for arrays in hits.iterate(list(branches), step_size=step_size, library="np"):
        mask = np.ones(len(arrays[branches[0]]), dtype=bool)
        for branch in branches:
            mask &= arrays[branch].astype(bool)
        selected += int(np.count_nonzero(mask))
    return selected


def summarize(path: Path, branches: tuple[str, ...], step_size: str) -> tuple[int, float, int]:
    with uproot.open(path) as root_file:
        if "RunMetadata" not in root_file:
            raise ValueError("missing RunMetadata tree; rerun with the metadata-enabled executable")
        if "MuonHits" not in root_file:
            raise ValueError("missing MuonHits tree")
        jobs, live_time_s = metadata_live_time(root_file["RunMetadata"])
        selected = count_coincidences(root_file["MuonHits"], branches, step_size)
    return jobs, live_time_s, selected


def print_result(label: str, jobs: int, live_time_s: float, selected: int,
                 target_events: int) -> None:
    days = live_time_s / 86400.0
    rate_hz = selected / live_time_s if live_time_s > 0.0 else 0.0
    rate_day = rate_hz * 86400.0
    target_days = target_events / rate_day if rate_day > 0.0 else math.inf
    poisson_percent = 100.0 / math.sqrt(selected) if selected > 0 else math.inf
    print(f"\n{label}")
    print(f"  merged jobs:              {jobs}")
    print(f"  selected events:           {selected}")
    print(f"  equivalent live time:      {days:.6f} days ({live_time_s:.3f} s)")
    print(f"  selected rate:             {rate_day:.8f} events/day ({rate_hz:.12g} Hz)")
    print(f"  days for {target_events:,} events:    {target_days:.6f}")
    print(f"  Poisson uncertainty now:   {poisson_percent:.3f}%")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("root_files", nargs="+", type=Path,
                        help="ROOT file(s), including hadd-merged files")
    parser.add_argument("--selection", choices=SELECTIONS, default="in2-out123",
                        help="built-in GEM coincidence (default: in2-out123)")
    parser.add_argument("--branches", nargs="+", metavar="BRANCH",
                        help="override --selection with GEM-valid branch names")
    parser.add_argument("--target-events", type=int, default=10_000,
                        help="report days required for this many selected events (default: 10000)")
    parser.add_argument("--step-size", default="50 MB",
                        help="uproot chunk size while counting hits (default: '50 MB')")
    args = parser.parse_args()

    if args.target_events <= 0:
        parser.error("--target-events must be positive")
    if args.branches:
        selection_label = ", ".join(args.branches)
        branches = tuple(args.branches)
    else:
        selection_label, branches = SELECTIONS[args.selection]

    print("GEM selection: " + selection_label)
    total_jobs = 0
    total_live_time_s = 0.0
    total_selected = 0
    for path in args.root_files:
        try:
            jobs, live_time_s, selected = summarize(path, branches, args.step_size)
        except (OSError, ValueError, KeyError) as error:
            print(f"error: {path}: {error}", file=sys.stderr)
            return 2
        print_result(str(path), jobs, live_time_s, selected, args.target_events)
        total_jobs += jobs
        total_live_time_s += live_time_s
        total_selected += selected

    if len(args.root_files) > 1:
        print_result("COMBINED", total_jobs, total_live_time_s,
                     total_selected, args.target_events)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
