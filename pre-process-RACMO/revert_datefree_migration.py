#!/usr/bin/env python3
"""
Reverse the date-free forcing migration done by migrate_existing_to_datefree.py.

Reads the manifest (datefree_migration_manifest.json) and:
  1. strips the stamped global attributes (dtobs from timeseries; start_ave_year /
     end_ave_year from averages), and
  2. renames every file back to its original date-bearing name.

Stripping dtobs from the multi-GB classic-netCDF timeseries files is a full
rewrite each, so it is parallelised across --workers (run via SLURM for big
domains). Stripping is idempotent: attributes already absent are skipped.

Dry-run by default. On a successful full revert the manifest is renamed to
<name>.reverted so it isn't applied twice.

    python3 revert_datefree_migration.py /path/to/<DOMAIN>_<FORCING>/input
    python3 revert_datefree_migration.py /path/to/.../input --apply --workers 32
"""

import argparse
import json
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

MANIFEST_NAME = "datefree_migration_manifest.json"
TS_ATTRS = ["dtobs"]
AVE_ATTRS = ["start_ave_year", "end_ave_year"]


def _has_global_attr(path, name):
    r = subprocess.run(["ncdump", "-h", str(path)], capture_output=True, text=True)
    return r.returncode == 0 and f":{name} = " in r.stdout


def _strip(path, names):
    """Delete the named global attrs that are present; returns True if it wrote."""
    present = [n for n in names if _has_global_attr(path, n)]
    if not present:
        return False
    cmd = ["ncatted", "-h", "-O"]
    for n in present:
        cmd += ["-a", f"{n},global,d,,"]
    cmd.append(str(path))
    subprocess.run(cmd, check=True, capture_output=True)
    return True


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_dir", help="Path to the <DOMAIN>_<FORCING>/input directory")
    p.add_argument("--apply", action="store_true", help="Actually revert (default: dry-run)")
    p.add_argument("--workers", type=int, default=8, help="Parallel attr-stripping workers (default 8)")
    args = p.parse_args()

    input_dir = Path(args.input_dir)
    manifest = input_dir / MANIFEST_NAME
    if not manifest.exists():
        sys.exit(f"ERROR: manifest not found: {manifest}")
    if args.apply:
        for tool in ("ncatted", "ncdump"):
            if shutil.which(tool) is None:
                sys.exit(f"ERROR: {tool} not found on PATH. Run `module load nco` first.")

    records = json.loads(manifest.read_text())
    ts = [r for r in records if r["subdir"] == "timeseries"]
    ave = [r for r in records if r["subdir"] == "averages"]
    mode = "APPLY" if args.apply else "DRY-RUN"
    print(f"[{mode}] reverting {len(records)} record(s) under {input_dir}\n")

    if not args.apply:
        print(f"  would strip {TS_ATTRS} from {len(ts)} timeseries files,")
        print(f"  strip {AVE_ATTRS} from {len(ave)} averages files,")
        print(f"  and rename all {len(records)} files back to date-bearing names.")
        print("\nDry-run only. Re-run with --apply --workers N.")
        sys.exit(0)

    # 1) strip attrs while names are still date-free
    print(f"  stripping {TS_ATTRS} from timeseries ({args.workers} workers)...")
    ts_paths = [input_dir / r["subdir"] / r["new"] for r in ts]
    done = wrote = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_strip, pth, TS_ATTRS): pth for pth in ts_paths if pth.exists()}
        for fut in as_completed(futs):
            done += 1
            if fut.result():
                wrote += 1
            if done % 50 == 0 or done == len(futs):
                print(f"    {done}/{len(futs)} processed, {wrote} rewritten", flush=True)
    for r in ave:
        pth = input_dir / r["subdir"] / r["new"]
        if pth.exists():
            _strip(pth, AVE_ATTRS)

    # 2) rename back to the original date-bearing names
    renamed = skipped = 0
    for r in records:
        d = input_dir / r["subdir"]
        new_path, old_path = d / r["new"], d / r["old"]
        if new_path.exists():
            if old_path.exists():
                sys.exit(f"ERROR: original already exists, aborting: {old_path}")
            new_path.rename(old_path)
            renamed += 1
        elif old_path.exists():
            skipped += 1  # already reverted

    print(f"\n  renamed back {renamed}, {skipped} already reverted")
    done_manifest = manifest.with_name(manifest.name + ".reverted")
    manifest.rename(done_manifest)
    print(f"  Manifest consumed -> {done_manifest.name}")
