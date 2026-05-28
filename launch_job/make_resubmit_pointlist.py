#!/usr/bin/env python3
"""
Check which points are fully complete (1D + 2D + 2Ddetail output files all exist)
and write a new pointlist containing only the remaining points.
Prints the count of remaining points to stdout — used by submit_job.sh.
"""
import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir',      required=True, help='Path to output/ directory')
    parser.add_argument('--pointlist',       required=True, help='Current pointlist file')
    parser.add_argument('--next-pointlist',  required=True, help='Path to write remaining points')
    parser.add_argument('--prefix',          required=True, help='Output file prefix, e.g. FGRN055_era055')
    args = parser.parse_args()

    with open(args.pointlist) as f:
        points = [int(line.strip()) for line in f if line.strip()]

    remaining = []
    for pt in points:
        f1d   = os.path.join(args.output_dir, f"{args.prefix}_1D_{pt}.nc")
        f2d   = os.path.join(args.output_dir, f"{args.prefix}_2D_{pt}.nc")
        f2ddet = os.path.join(args.output_dir, f"{args.prefix}_2Ddetail_{pt}.nc")
        if not (os.path.exists(f1d) and os.path.exists(f2d) and os.path.exists(f2ddet)):
            remaining.append(pt)

    with open(args.next_pointlist, 'w') as f:
        for pt in remaining:
            f.write(f"{pt}\n")

    print(len(remaining))


if __name__ == '__main__':
    main()
