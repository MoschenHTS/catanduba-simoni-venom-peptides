#!/usr/bin/env python3
"""
Rewrite the absolute #include paths inside GROMACS topologies to relative ones.

WHY THIS IS A SEPARATE SCRIPT, NOT DONE IN PLACE
The committed md/<system>/{topol,system}.top files are PRIMARY DATA - they are the exact
topologies the published 200 ns x 8 replica simulations were run with. Silently editing
them would make the archive no longer match what was actually simulated. So they ship
verbatim, absolute paths and all, and this script converts a COPY on demand.

WHAT'S ABSOLUTE IN THEM, AND WHY
  - forcefield.itp, tip3p.itp, ions.itp : written by pdb2gmx, which records the full path
    to whichever force field it resolved (here md/charmm36-jul2022.ff/)
  - POPC.itp : written by 03_embed.py (fixed 2026-08 to emit a relative path for new runs;
    the committed topologies predate that fix)

THE FAILURE MODE THIS PREVENTS
grompp does not error on an absolute include that happens to exist. If a stale copy of
the original tree is still present on the machine, grompp will silently read the force
field from THERE while every other script reads the new location. If the two trees ever
diverge, the simulation is quietly wrong with no warning. Convert before reusing.

Usage:
    python 13_normalize_top_paths.py <file.top> [...]        # rewrite in place
    python 13_normalize_top_paths.py --check <file.top> [...] # report only, exit 1 if any found
"""
import os
import re
import sys

OLD_BASE = "/mnt/SSD3/henrique/c_simoni/boltz_labpc/md"

# include-path -> path relative to the .top file's own directory (md/<system>/)
RELATIVE = {
    "charmm36-jul2022.ff/forcefield.itp": "../charmm36-jul2022.ff/forcefield.itp",
    "charmm36-jul2022.ff/tip3p.itp": "../charmm36-jul2022.ff/tip3p.itp",
    "charmm36-jul2022.ff/ions.itp": "../charmm36-jul2022.ff/ions.itp",
    "bilayer/POPC.itp": "../bilayer/POPC.itp",
}

INCLUDE = re.compile(r'#include\s+"([^"]+)"')


def convert(path, check_only=False):
    text = open(path, encoding="utf-8").read()
    hits, out = [], text
    for m in INCLUDE.finditer(text):
        inc = m.group(1)
        if not os.path.isabs(inc):
            continue
        hits.append(inc)
        if check_only:
            continue
        tail = inc[len(OLD_BASE):].lstrip("/") if inc.startswith(OLD_BASE) else None
        rel = RELATIVE.get(tail)
        if rel is None:
            # Unknown absolute include: fall back to a path relative to the .top file,
            # rather than guessing - and say so loudly.
            rel = os.path.relpath(inc, os.path.dirname(os.path.abspath(path)))
            print(f"  ! {path}: unrecognised absolute include {inc}\n"
                  f"    -> using computed relative path {rel} (VERIFY THIS)")
        out = out.replace(f'#include "{inc}"', f'#include "{rel}"')
    if check_only:
        for h in hits:
            print(f"  {path}: {h}")
        return len(hits)
    if hits:
        open(path, "w", encoding="utf-8").write(out)
        print(f"  {path}: rewrote {len(hits)} absolute include(s)")
    else:
        print(f"  {path}: already relative, unchanged")
    return len(hits)


if __name__ == "__main__":
    args = sys.argv[1:]
    check = "--check" in args
    files = [a for a in args if a != "--check"]
    if not files:
        raise SystemExit(__doc__)
    total = sum(convert(f, check_only=check) for f in files)
    if check:
        print(f"\n{total} absolute include(s) found across {len(files)} file(s)")
        sys.exit(1 if total else 0)
