#!/usr/bin/env python3
"""
Migrate existing IMAU-FDM forcing files to the date-free naming AND stamp the
same global metadata that fresh preprocessing writes, so migrated files are
equivalent to freshly-produced ones.

Timeseries:  <var>_<domain>_<forcing>_<YYYY>-<YYYY>_p<N>.nc
          -> <var>_<domain>_<forcing>_p<N>.nc          + global attr dtobs (EVERY file)

Averages:    <var>_<domain>_<forcing>-<YYYY>_<YYYY>-<YYYY>_ave.nc
          -> <var>_<domain>_<forcing>_ave.nc           + global attrs start_ave_year,
                                                          end_ave_year (parsed from old name)

Renames are instant and done serially. Stamping dtobs on the (multi-GB, classic
netCDF) timeseries files requires a full file rewrite each, so it is parallelised
across --workers and is best run as a SLURM job (see migrate_datefree.sbatch).
Stamping is idempotent: files that already carry the attribute are skipped, so the
job is safe to re-run.

On --apply a manifest (datefree_migration_manifest.json) recording every rename is
written into the input dir, so the change can be undone with
revert_datefree_migration.py.

    # dry-run (default)
    python3 migrate_existing_to_datefree.py /path/to/<DOMAIN>_<FORCING>/input

    # apply with 32 parallel stamping workers (needs `module load nco`)
    python3 migrate_existing_to_datefree.py /path/to/.../input --apply --workers 32
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Forcing timestep [s]; matches DTOBS_SECONDS in the preprocessing config (3-hourly ERA5).
DTOBS_SECONDS = 10800
MANIFEST_NAME = "datefree_migration_manifest.json"

# <base>_<YYYY>-<YYYY>_p<N>.nc       e.g. evap_FGRN055_era5_1939-2025_p10.nc
TS_RE = re.compile(r"^(?P<base>.+)_(?P<start>\d{4})-(?P<end>\d{4})_p(?P<part>\d+)\.nc$")
# <base>-<YYYY>_<YYYY>-<YYYY>_ave.nc e.g. evap_FGRN055_era5-1939_1940-1970_ave.nc
AVE_RE = re.compile(r"^(?P<base>.+)-(?P<tsstart>\d{4})_(?P<avestart>\d{4})-(?P<aveend>\d{4})_ave\.nc$")


def _has_global_attr(path, name):
    r = subprocess.run(["ncdump", "-h", str(path)], capture_output=True, text=True)
    return r.returncode == 0 and f":{name} = " in r.stdout


def _ncatted_set(path, attrs):
    """Set integer global attrs (-h: don't touch the history attribute)."""
    cmd = ["ncatted", "-h", "-O"]
    for name, val in attrs:
        cmd += ["-a", f"{name},global,o,l,{val}"]
    cmd.append(str(path))
    subprocess.run(cmd, check=True, capture_output=True)


def rename_phase(input_dir):
    """Rename date-bearing files to date-free; stamp ave-year attrs on averages
    (small files, done here serially). Returns list of rename records."""
    ts_dir, ave_dir = input_dir / "timeseries", input_dir / "averages"
    records = []

    for f in sorted(ts_dir.glob("*_p*.nc")):
        m = TS_RE.match(f.name)
        if not m:
            continue  # already date-free
        new = ts_dir / f"{m['base']}_p{m['part']}.nc"
        if new.exists():
            sys.exit(f"ERROR: target already exists, aborting: {new}")
        f.rename(new)
        records.append({"subdir": "timeseries", "old": f.name, "new": new.name})

    for f in sorted(ave_dir.glob("*_ave.nc")):
        m = AVE_RE.match(f.name)
        if not m:
            continue  # already date-free
        new = ave_dir / f"{m['base']}_ave.nc"
        if new.exists():
            sys.exit(f"ERROR: target already exists, aborting: {new}")
        _ncatted_set(f, [("start_ave_year", int(m["avestart"])),
                         ("end_ave_year", int(m["aveend"]))])
        f.rename(new)
        records.append({"subdir": "averages", "old": f.name, "new": new.name})

    return records


def _stamp_dtobs(path):
    """Idempotently stamp dtobs; returns True if it had to write the file."""
    if _has_global_attr(path, "dtobs"):
        return False
    _ncatted_set(path, [("dtobs", DTOBS_SECONDS)])
    return True


def stamp_phase(input_dir, workers):
    """Stamp dtobs on every date-free timeseries file (parallel, idempotent)."""
    ts_dir = input_dir / "timeseries"
    files = [f for f in sorted(ts_dir.glob("*_p*.nc")) if not TS_RE.match(f.name)]
    total = len(files)
    print(f"  stamping dtobs on {total} timeseries file(s) with {workers} worker(s)...")
    done = stamped = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_stamp_dtobs, f): f for f in files}
        for fut in as_completed(futs):
            done += 1
            try:
                if fut.result():
                    stamped += 1
            except subprocess.CalledProcessError as e:
                sys.exit(f"ERROR stamping {futs[fut].name}: {e}")
            if done % 50 == 0 or done == total:
                print(f"    {done}/{total} processed, {stamped} newly stamped", flush=True)
    print(f"  dtobs present on all {total} timeseries files ({stamped} newly stamped)\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_dir", help="Path to the <DOMAIN>_<FORCING>/input directory")
    p.add_argument("--apply", action="store_true", help="Actually rename + stamp (default: dry-run)")
    p.add_argument("--workers", type=int, default=8, help="Parallel dtobs-stamping workers (default 8)")
    args = p.parse_args()

    input_dir = Path(args.input_dir)
    for d in (input_dir / "timeseries", input_dir / "averages"):
        if not d.is_dir():
            sys.exit(f"ERROR: not a directory: {d}")
    if args.apply:
        for tool in ("ncatted", "ncdump"):
            if shutil.which(tool) is None:
                sys.exit(f"ERROR: {tool} not found on PATH. Run `module load nco` first.")

    mode = "APPLY" if args.apply else "DRY-RUN"
    print(f"[{mode}] migrating forcing files under {input_dir}\n")

    if not args.apply:
        ts, ave = input_dir / "timeseries", input_dir / "averages"
        n_ren_ts = sum(1 for f in ts.glob("*_p*.nc") if TS_RE.match(f.name))
        n_ren_ave = sum(1 for f in ave.glob("*_ave.nc") if AVE_RE.match(f.name))
        n_df_ts = sum(1 for f in ts.glob("*_p*.nc") if not TS_RE.match(f.name))
        print(f"  timeseries to rename: {n_ren_ts}; averages to rename: {n_ren_ave}")
        print(f"  date-free timeseries present: {n_df_ts} (each checked for dtobs, stamped if absent)")
        print("\nDry-run only. Re-run with --apply --workers N.")
        sys.exit(0)

    records = rename_phase(input_dir)
    if records:
        manifest = input_dir / MANIFEST_NAME
        existing = json.loads(manifest.read_text()) if manifest.exists() else []
        existing.extend(records)
        manifest.write_text(json.dumps(existing, indent=2) + "\n")
        print(f"  renamed {len(records)} file(s); manifest -> {manifest}\n")
    else:
        print("  no date-bearing files to rename (already done)\n")

    stamp_phase(input_dir, args.workers)
    print("Migration complete.")
