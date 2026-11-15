#!/usr/bin/env python3
"""Check the ensemble of examples/fault_system_uncertainty.yaml (STX-TN-13, Prompt 41 definition
of done): N members recorded, each with a status and its checks; ranges state N and the seed; and
the first-order indices rank the parameters as the analytical case predicts:
  - Max throw F3 = F3's throw: faults.F3.throw.max first;
  - Max throw F2 = F2's throw (members with the splay): faults.F2.throw.max first;
  - Max throw F1 = F1's throw plus the splay's throw transferred onto F1 beyond the branch line:
    the splay scenario (with or without F2, variance 0.21 x 50^2) first, then F1's own throw
    (15^2), then F3's throw (through the 0.6 correlation); the well top not distinguishable.
Usage: check_ensemble.py <ensemble output folder> <members> <seed>
"""
import csv
import pathlib
import sys


def main(folder, members, seed):
    e = pathlib.Path(folder) / "ensemble"
    bad = []
    rows = list(csv.DictReader(open(e / "members.csv")))
    mem = [r for r in rows if r["member"] != "base"]
    if len(mem) != members:
        bad.append("%d members recorded, expected %d" % (len(mem), members))
    for r in mem:
        if r["status"] not in ("ok", "refused", "build_failed", "checks_failed"):
            bad.append("member %s has no status" % r["member"])
        if r["status"] != "ok" and not r["reason"]:
            bad.append("member %s failed without a reason" % r["member"])
    for r in csv.DictReader(open(e / "statistics.csv")):
        if r["members"] != str(members) or r["seed"] != str(seed):
            bad.append("statistics row %s does not state N and the seed" % r["key"])
    sob = {}
    for r in csv.DictReader(open(e / "sensitivity_sobol.csv")):
        sob.setdefault(r["key"], []).append((float(r["first_order_index"]), r["parameter"], r["distinguishable_from_0"]))
    for k in sob:
        sob[k].sort(reverse=True)
        print("%s: %s" % (k, ", ".join("%s %.3f" % (p, v) for v, p, _ in sob[k])))
    def order(key, expected):
        got = [p for _, p, _ in sob.get(key, [])][:len(expected)]
        if got != expected:
            bad.append("%s ranked %s, expected %s" % (key, got, expected))
    order("Max throw F3", ["faults.F3.throw.max"])
    order("Max throw F2", ["faults.F2.throw.max"])
    order("Max throw F1", ["scenario splay", "faults.F1.throw.max", "faults.F3.throw.max"])
    for v, p, d in sob.get("Max throw F1", []):
        if p == "wells.W1.tops.H2" and d == "yes":
            bad.append("the well top should not be distinguishable from 0 for Max throw F1 (%.3f)" % v)
    ok = sum(r["status"] == "ok" for r in mem)
    print("%d members, %d ok, seed %d" % (len(mem), ok, seed))
    if bad:
        print("\n".join(bad))
        return 1
    print("ensemble checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3])))
