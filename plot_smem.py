#!/usr/bin/env python3
"""
Haplotype x site matrix with SMEMs colored in place.
  red    = cell covered only by reciprocal SMEMs (i->j and j->i both reported)
  blue   = cell covered only by one-sided SMEMs
  purple = covered by both kinds
  white/gray = not in any SMEM (gray = allele 1)
--paint both  : each SMEM colors mh and neighbor rows (default)
--paint query : each SMEM colors only the mh row (the haplotype that reported it)

Usage:
  python plot_smem_matrix.py Algo4/s100i50.rnd.vcf algo4.out --min-len 20
  python plot_smem_matrix.py data.vcf algo4.out --min-len 20 --region 0:100 --haps 0-39 --out m.png
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


def read_vcf(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        r = []
        for g in line.split()[9:]:
            a, b = g.split(":")[0].split("|")
            r += [int(a), int(b)]
        rows.append(r)
    return np.array(rows, dtype=np.int8).T  # haplotypes x sites


def read_smems(path, min_len):
    out = []
    for line in open(path):
        p = list(map(int, line.split()))
        if len(p) == 4:                      # current C++ output: mh nb k+1 len
            i, j, E, L = p
            s, e = E - 1 - L, E - 2          # inclusive 0-based interval
        elif len(p) >= 5:                    # mh nb start end len
            i, j, s, e, L = p[:5]
        else:
            continue
        if L >= min_len and i != j:
            out.append((i, j, s, e))
    keys = set(out)
    return [(i, j, s, e, (j, i, s, e) in keys) for i, j, s, e in out]


def parse_haps(txt, M):
    if not txt:
        return list(range(M))
    idx = []
    for part in txt.split(","):
        a, _, b = part.partition("-")
        idx += range(int(a), int(b or a) + 1)
    return [h for h in idx if h < M]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf"); ap.add_argument("smems")
    ap.add_argument("--min-len", type=int, default=1)
    ap.add_argument("--region", help="sites a:b (0-based, b exclusive)")
    ap.add_argument("--haps", help="e.g. 0-39 or 0,5,10-20")
    ap.add_argument("--paint", choices=["both", "query"], default="both")
    ap.add_argument("--out")
    args = ap.parse_args()

    X = read_vcf(args.vcf)
    M, N = X.shape
    a, b = map(int, args.region.split(":")) if args.region else (0, N)

    R = np.zeros((M, N), bool)   # covered by a reciprocal SMEM
    B = np.zeros((M, N), bool)   # covered by a one-sided SMEM
    smems = read_smems(args.smems, args.min_len)
    for i, j, s, e, rec in smems:
        for h in ((i, j) if args.paint == "both" else (i,)):
            (R if rec else B)[h, s:e + 1] = True
    n_one = sum(not t[4] for t in smems)
    print(f"SMEMs (len >= {args.min_len}): {len(smems)} | reciprocal {len(smems) - n_one} "
          f"| one-sided {n_one} | blue cells {B.sum()} (overlapping red: {(B & R).sum()})")
    # 0 none, 1 one-sided only, 2 reciprocal only, 3 both
    C = (B & ~R) * 1 + (R & ~B) * 2 + (R & B) * 3

    haps = parse_haps(args.haps, M)
    sub_c = C[haps, a:b]
    sub_x = X[haps, a:b]
    # color code: none&0, none&1, blue, red, purple
    img = np.where(sub_c == 0, sub_x, sub_c + 1)
    cmap = ListedColormap(["white", "#d9d9d9", "#5b8fd9", "#e05a5a", "#a46fd1"])

    nh, ns = sub_x.shape
    fig, ax = plt.subplots(figsize=(min(24, 2 + ns * 0.25), min(30, 1.5 + nh * 0.25)))
    ax.imshow(img, cmap=cmap, vmin=0, vmax=4, aspect="auto", interpolation="nearest",
              extent=(a - .5, b - .5, nh - .5, -.5))
    if nh * ns <= 6000:                       # write the 0/1 digits when readable
        for r in range(nh):
            for c in range(ns):
                ax.text(a + c, r, sub_x[r, c], ha="center", va="center", fontsize=7)
        ax.set_xticks(np.arange(a - .5, b, 1), minor=True)
        ax.set_yticks(np.arange(-.5, nh, 1), minor=True)
        ax.grid(which="minor", color="white", lw=0.5)
        ax.tick_params(which="minor", length=0)
    ax.set_yticks(range(nh)); ax.set_yticklabels(haps, fontsize=7)
    ax.set_xlabel("site"); ax.set_ylabel("haplotype")
    ax.set_title(f"SMEMs on haplotype x site matrix (length >= {args.min_len}, paint={args.paint})")
    ax.legend(handles=[Patch(color="#e05a5a", label="reciprocal SMEM"),
                       Patch(color="#5b8fd9", label="one-sided SMEM"),
                       Patch(color="#a46fd1", label="both"),
                       Patch(facecolor="white", edgecolor="gray", label="no SMEM")],
              loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, frameon=False)
    plt.tight_layout()
    if args.out:
        plt.savefig(args.out, dpi=200, bbox_inches="tight"); print("saved", args.out)
    else:
        plt.show()


if __name__ == "__main__":
    main()