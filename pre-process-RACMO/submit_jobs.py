#!/usr/bin/env python3
"""
Submit SLURM pre-processing jobs for RACMO forcing files.
Replaces run_sc_job.sc + all_make*.sc

One sbatch job is submitted per variable per step. When multiple steps are
requested, each step waits for all jobs from the previous step to complete
successfully (SLURM afterok dependency), so years → timeseries → averages
run in the correct order without manual intervention.

Run from within the domain subdirectory (which contains config.py):

    cd FGRN055_era055
    python3 ../submit_jobs.py all

Usage examples:
    # Full pipeline for all variables
    python3 ../submit_jobs.py all

    # Only create yearly files
    python3 ../submit_jobs.py years

    # Create timeseries without re-creating parts (e.g. extending timeseries)
    python3 submit_jobs.py timeseries --no-parts

    # years + timeseries, skip averages (e.g. when extending an existing run)
    python3 submit_jobs.py years timeseries --skip-averages

    # Extend timeseries with new years — create parts only (safe, non-destructive)
    python3 submit_jobs.py years timeseries --extend-start 2024 --extend-end 2025 --extend-mode parts-only

    # Extend timeseries with new years — append to existing timeseries files
    python3 submit_jobs.py years timeseries --extend-start 2024 --extend-end 2025 --extend-mode append

    # Subset of variables
    python3 submit_jobs.py all --vars precip snowfall

    # Check what would be submitted without actually submitting
    python3 submit_jobs.py all --dry-run
"""

import argparse
import re
import subprocess
from pathlib import Path

import config

SCRIPT_DIR = Path(__file__).parent   # pre-process-RACMO/ — where scripts live
DOMAIN_DIR = Path.cwd()              # domain subdir the user ran from (contains config.py)

VALID_JOBS = ["years", "timeseries", "averages"]


def _build_timeseries_args(args):
    """Build extra CLI args string for make_fdm_timeseries.py from parsed args."""
    parts = []
    if args.no_parts:
        parts.append("--no-parts")
    if args.extend_start:
        parts += ["--extend-start", str(args.extend_start),
                  "--extend-end",   str(args.extend_end),
                  "--extend-mode",  args.extend_mode]
    return " ".join(parts)


def create_job_script(jobname, varname, extra_args=""):
    """Write an sbatch job file and return its path."""
    script_map = {
        "years":      ("make_fdm_years.py",      varname),
        "timeseries": ("make_fdm_timeseries.py",  f"{varname} {extra_args}".strip()),
        "averages":   ("make_fdm_averages.py",    varname),
    }
    py_script, call_args = script_map[jobname]

    config.JOBFILE_DIR.mkdir(parents=True, exist_ok=True)
    config.LOGFILE_DIR.mkdir(parents=True, exist_ok=True)

    jobfile = config.JOBFILE_DIR / f"preprocess-RACMO_{jobname}_{varname}"
    logfile = config.LOGFILE_DIR / f"{jobname}_{varname}_preprocess.log"

    script = f"""#!/bin/bash

#SBATCH -q {config.SLURM_QUEUE}
#SBATCH -J {jobname}_{varname}_preprocess
#SBATCH --time={config.SLURM_TIME}
#SBATCH -o {logfile}
#SBATCH --mem-per-cpu={config.SLURM_MEM}

module load nco

echo "Job: {jobname}  Variable: {varname}"
echo "TS years:  {config.TS_START_YEAR}-{config.TS_END_YEAR}"
echo "Ave years: {config.AVE_START_YEAR}-{config.AVE_END_YEAR}"
echo "Lon bands: {config.NUM_LONG_BANDS} x {config.CELL_WIDTH}"

cd {DOMAIN_DIR}
python3 {SCRIPT_DIR}/{py_script} {call_args}
"""
    jobfile.write_text(script)
    return jobfile


def _parse_job_id(sbatch_output):
    """Extract integer job ID from 'Submitted batch job 12345'."""
    m = re.search(r"(\d+)", sbatch_output)
    return m.group(1) if m else None


def submit(jobs, varnames, dry_run=False, extra_args=None):
    """
    Submit jobs step by step. Each step collects the job IDs from the previous
    step and passes them as afterok dependencies, so the order years →
    timeseries → averages is enforced automatically.
    """
    extra_args = extra_args or {}
    prev_job_ids = []  # job IDs from the previous step

    for jobname in jobs:
        current_job_ids = []
        dep_flag = (
            f"--dependency=afterok:{':'.join(prev_job_ids)}"
            if prev_job_ids else None
        )

        for varname in varnames:
            jobfile = create_job_script(jobname, varname, extra_args.get(jobname, ""))

            if dry_run:
                dep_str = f" (after: {','.join(prev_job_ids)})" if prev_job_ids else ""
                print(f"[dry-run] {jobname}/{varname}{dep_str}: {jobfile}")
                continue

            cmd = ["sbatch"]
            if dep_flag:
                cmd.append(dep_flag)
            cmd.append(str(jobfile))

            result = subprocess.run(cmd, capture_output=True, text=True, cwd=config.JOBFILE_DIR)
            output = result.stdout.strip() or result.stderr.strip()
            job_id = _parse_job_id(output)
            if job_id:
                current_job_ids.append(job_id)
            print(f"Submitted {jobname}/{varname}: {output}")

        if not dry_run:
            prev_job_ids = current_job_ids


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Submit RACMO pre-processing jobs to SLURM.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "jobs",
        nargs="+",
        choices=VALID_JOBS + ["all"],
        metavar="{years,timeseries,averages,all}",
        help="Which steps to submit ('all' expands to all three steps)",
    )
    parser.add_argument(
        "--vars",
        nargs="+",
        default=config.VARS,
        metavar="VAR",
        help=f"Variables to process (default: all {len(config.VARS)})",
    )
    parser.add_argument(
        "--no-parts",
        action="store_true",
        help="Pass --no-parts to timeseries jobs (skip re-creating lon-band parts)",
    )
    parser.add_argument(
        "--skip-averages",
        action="store_true",
        help="Omit the averages step (e.g. when extending an existing timeseries)",
    )
    parser.add_argument(
        "--extend-start",
        type=int,
        metavar="YEAR",
        help="First year of new RACMO data (triggers extend mode in timeseries jobs)",
    )
    parser.add_argument(
        "--extend-end",
        type=int,
        metavar="YEAR",
        help="Last year of new RACMO data",
    )
    parser.add_argument(
        "--extend-mode",
        choices=["parts-only", "append"],
        default="parts-only",
        help=(
            "parts-only: write new parts to parts_extend/ only; "
            "append: also append new parts to existing timeseries files "
            "(default: parts-only)"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write job files and print paths without submitting",
    )
    args = parser.parse_args()

    if (args.extend_start is None) != (args.extend_end is None):
        parser.error("--extend-start and --extend-end must both be provided together")

    jobs = []
    for j in args.jobs:
        if j == "all":
            jobs.extend(VALID_JOBS)
        else:
            jobs.append(j)
    # Deduplicate while preserving order
    seen = set()
    jobs = [j for j in jobs if not (j in seen or seen.add(j))]

    if args.skip_averages and "averages" in jobs:
        jobs.remove("averages")

    extra_args = {
        "timeseries": _build_timeseries_args(args),
    }

    submit(jobs, args.vars, dry_run=args.dry_run, extra_args=extra_args)
