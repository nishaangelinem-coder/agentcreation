#!/usr/bin/env python3
"""Merge per-worker synthesis results into results/syn_results.json."""
import json, glob, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
out = json.load(open("results/syn_results.json")) if os.path.exists("results/syn_results.json") else {}
for f in sorted(glob.glob("results/syn_[b-z].json")):
    out.update(json.load(open(f)))
json.dump(out, open("results/syn_results.json", "w"), indent=1)
print(len(out), "records")
