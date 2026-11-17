"""Draw the README figures of certified-bowen-york-spin from the data the repository ships.

Usage:  python figures/make_figures.py <repo root> [--png]
Needs Python 3.11 and matplotlib 3.11 (nothing else). Writes figures/*.svg next to this script
(or into ./figures when the script does not sit in a folder named figures). --png also writes PNG previews.

Every plotted number is read from shipped files (exact rationals "p/q" where the file stores them):
  results/p2e/result_p2e.json     22 certified nodes j = 2^k (chi_lo, chi_hi), the 23 certified upper-bound pieces
                                  (j_lo, j_hi, value), and the supremum window (X3, c)          [sealed test P2e]
  results/p2d/result_p2d.json     the earlier supremum window (X3_dec_down, c_dec_up)          [sealed test P2d]
  results/b3t/result_b3t.json     the extreme (mu = 0) throat: chi(0) window (chi_lo, chi_hi)  [sealed test B3T]
chi = J / M_ADM^2 and j = 4J/mu^2 (README, "How it works", Step 2).
"""
import json
import os
import sys
from fractions import Fraction as F

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

# ---------------------------------------------------------------- brand style (figures/STYLE.md)
BG, BORDER, GRID, AXC, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
CAT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]
W_IN = 12.5  # 900 pt viewBox width

plt.rcParams.update({
    "svg.fonttype": "none", "svg.hashsalt": "moki-julio",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 16, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "axes.facecolor": BG, "figure.facecolor": "none", "savefig.facecolor": "none",
    "axes.edgecolor": AXC, "axes.labelcolor": TXT, "xtick.color": AXC, "ytick.color": AXC,
    "text.color": TXT, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.fontsize": 14, "lines.linewidth": 2.0, "axes.axisbelow": True,
})
PNG = "--png" in sys.argv
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def new_fig(h_in):
    fig = plt.figure(figsize=(W_IN, h_in))
    fig.patches.append(FancyBboxPatch((0.004, 0.006), 0.992, 0.988, boxstyle="round,pad=0,rounding_size=0.015",
                                      transform=fig.transFigure, facecolor=BG, edgecolor=BORDER, linewidth=1.5,
                                      zorder=-10))
    return fig


def title(fig, main, sub):
    h = fig.get_figheight()
    fig.text(0.035, 1 - 0.42 / h, main, fontsize=21, fontweight="bold", color=TXT, va="top")
    fig.text(0.035, 1 - 0.88 / h, sub, fontsize=15, color=AXC, va="top")


def save(fig, out, name):
    path = os.path.join(out, name)
    fig.savefig(path, format="svg", metadata={"Date": None, "Creator": None})
    with open(path, encoding="utf-8") as fh:  # mathtext pieces name one font only; give them the full stack
        svg = fh.read().replace("font-family: 'Segoe UI';", "font-family: 'Segoe UI', 'Helvetica', 'Arial', sans-serif;")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(svg)
    if PNG:
        prev = os.path.join(out, os.pardir, "_preview")
        os.makedirs(prev, exist_ok=True)
        fig.savefig(os.path.join(prev, name.replace(".svg", ".png")), format="png", dpi=144)
    plt.close(fig)
    print("wrote figures/" + name)


def load(root, rel):
    with open(os.path.join(root, rel), encoding="utf-8") as fh:
        return json.load(fh)


def sci(x, nd=1):
    m, e = ("%.*e" % (nd, float(x))).split("e")
    return "%s × 10%s" % (m, str(int(e)).translate(SUP))


def dec(x, nd, direction):
    """exact rational -> decimal string with nd decimals, rounded outward ('down' for lower ends, 'up' for upper)."""
    q = x * 10 ** nd
    n = q.numerator // q.denominator if direction == "down" else -((-q.numerator) // q.denominator)
    s = str(n).rjust(nd + 1, "0")
    return s[:-nd] + "." + s[-nd:]


def pow2(j):
    k = F(j).numerator.bit_length() - 1
    assert F(j) == 2 ** k
    return "2" + str(k).translate(SUP) if k > 1 else str(2 ** k)


def data(root):
    p2e = load(root, "results/p2e/result_p2e.json")
    nodes = [(F(n["j"]), F(n["chi_lo"]), F(n["chi_hi"])) for n in p2e["nodes"]]
    assert all(n["point_certified"] for n in p2e["nodes"]) and len(nodes) == 22
    pieces = []
    for p in p2e["pieces"]:
        hi = None if p["j_hi"] == "oo" else F(p["j_hi"])
        pieces.append((F(p["j_lo"]), hi, F(p["value"]), p["label"]))
    win = load(root, "results/b3t/result_b3t.json")["window"]
    chi0 = (F(win["chi_lo"]), F(win["chi_hi"]))
    sup_e = (F(p2e["X3"]), F(p2e["c"]))
    p2d = load(root, "results/p2d/result_p2d.json")
    sup_d = (F(p2d["X3"]), F(p2d["c"]))
    return nodes, pieces, chi0, sup_e, sup_d


# ---------------------------------------------------------------- figure 1: the climb to the ceiling (hero)
def fig_ceiling(root, out):
    nodes, pieces, chi0, sup_e, sup_d = data(root)
    fig = new_fig(9.4)
    title(fig, "A spin ceiling, proved",
          "Spin J/M² of the Bowen-York black hole as the bare mass shrinks (j grows); every value certified")

    ax = fig.add_axes([0.10, 0.475, 0.86, 0.36])
    js = [float(j) for j, _, _ in nodes]
    # certified upper bound on each stretch of j (one horizontal step per piece; the tail runs to infinity)
    jmax = js[-1] * 3.2
    for lo, hi, v, lab in pieces:
        if lo == 0:
            continue
        x1 = float(hi) if hi is not None else jmax
        ax.plot([float(lo), x1], [float(v), float(v)], color=CAT[0], lw=2.0, solid_capstyle="butt", zorder=2)
    for j, lo, hi in nodes:
        ax.plot([float(j), float(j)], [float(lo), float(hi)], color=ACC, lw=2.0, zorder=3)
        ax.plot(float(j), float((lo + hi) / 2), "o", ms=8, color=ACC, zorder=4)
    c0 = float(sum(chi0) / 2)
    ax.axhline(1.0, color=AXC, lw=1, ls=(0, (5, 4)))
    ax.text(1.05, 1.008, "1: the most a settled black hole can spin (for comparison)", color=AXC, fontsize=13,
            va="bottom")
    ax.axhline(c0, color=CAT[1], lw=1.6, zorder=1)
    ax.annotate("the ceiling: spin of the endless throat (μ = 0)", (2 ** 9, c0), xytext=(2 ** 7.2, 0.80),
                color=CAT[1], fontsize=14, va="center", arrowprops=dict(arrowstyle="-", color=CAT[1], lw=1.2))
    ax.set_xscale("log", base=2)
    ax.set_xlim(0.8, jmax)
    ax.set_ylim(0.2, 1.09)
    ax.set_ylabel("spin J/M²")
    ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=2, numticks=12))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: (str(int(v)) if v < 64 else "2" + str(int(round(__import__("math").log2(v)))).translate(SUP))))
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.legend(handles=[Line2D([], [], color=ACC, marker="o", ms=8, label="certified value at j = 1, 2, 4, …, 2²¹"),
                       Line2D([], [], color=CAT[0], lw=2, label="certified upper bound for every j in each stretch")],
              loc="lower right", fontsize=14)

    bx = fig.add_axes([0.10, 0.09, 0.86, 0.32], sharex=ax)
    gaps = []
    for j, lo, hi in nodes:
        g_lo, g_hi = chi0[0] - hi, chi0[1] - lo
        assert g_lo > 0
        gaps.append((j, g_lo, g_hi))
        bx.plot([float(j), float(j)], [float(g_lo), float(g_hi)], color=ACC, lw=2.0, zorder=3)
        bx.plot(float(j), float((g_lo + g_hi) / 2), "o", ms=8, color=ACC, zorder=4)
    bx.set_yscale("log")
    bx.set_ylabel("distance below the\nceiling (log scale)")
    bx.set_xlabel("j = 4J/μ²  (log scale; large j = small bare mass μ)")
    last = gaps[-1]
    bx.annotate("at j = 2²¹: at least %s below\n(the ceiling is never reached)" % sci(last[1]),
                (float(last[0]), float(last[1])), xytext=(2 ** 0.5, 1.5e-8), color=TXT, fontsize=14, va="center",
                arrowprops=dict(arrowstyle="-", color=AXC, lw=1.2))
    save(fig, out, "spin-ceiling.svg")
    return dict(n_nodes=len(nodes), first=(str(nodes[0][0]), float(nodes[0][1]), float(nodes[0][2])),
                last=(str(nodes[-1][0]), float(nodes[-1][1]), float(nodes[-1][2])),
                chi0=[float(x) for x in chi0], gap_first=[float(gaps[0][1]), float(gaps[0][2])],
                gap_last=[float(last[1]), float(last[2])], n_pieces=len(pieces),
                max_node_width=max(float(h - l) for _, l, h in nodes))


# ---------------------------------------------------------------- figure 2: pinning the ceiling
def fig_windows(root, out):
    nodes, pieces, chi0, sup_e, sup_d = data(root)
    ref = F(928200, 10 ** 6)
    unit = F(1, 10 ** 9)
    rows = [("earlier certificate (P2d):\nthe ceiling lies in this window", sup_d, CAT[0]),
            ("final certificate (P2e):\nthe ceiling lies in this window", sup_e, CAT[0]),
            ("the endless throat itself (B3T):\nits spin lies in this window", chi0, CAT[1])]
    fig = new_fig(6.6)
    title(fig, "Pinning the ceiling",
          "Each certified window sits inside the one before; position shown as J/M² − 0.9282, in billionths")
    ax = fig.add_axes([0.37, 0.14, 0.59, 0.62])
    info = []
    for i, (lab, (lo, hi), col) in enumerate(rows):
        y = len(rows) - 1 - i
        x0, x1 = float((lo - ref) / unit), float((hi - ref) / unit)
        ax.add_patch(Rectangle((x0, y - 0.28), x1 - x0, 0.56, facecolor=col, alpha=0.30, edgecolor=col, lw=1))
        if x1 - x0 < 0.5:
            ax.plot([x0, x0], [y - 0.32, y + 0.32], color=col, lw=2.4)
        ax.text(-0.03, y, lab, fontsize=14, color=TXT, va="center", ha="right",
                transform=blended_transform_factory(ax.transAxes, ax.transData))
        nd = 14 if i == 2 else 10
        txt = "[%s, %s]" % (dec(lo, nd, "down"), dec(hi, nd, "up"))
        if i == 2:
            txt += "\nwidth %s" % sci(hi - lo)
        ax.text(x1 + 6, y, txt, color=col, fontsize=13, va="center")
        info.append(dict(row=lab.split("\n")[0], lo=str(lo.numerator * 10 ** 22 // lo.denominator),
                         hi=str(-((-hi.numerator * 10 ** 22) // hi.denominator)), width=float(hi - lo)))
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlim(float((sup_d[0] - ref) / unit) - 40, float((sup_d[1] - ref) / unit) + 330)
    ax.set_xlabel("J/M² − 0.9282  (units of 10⁻⁹)")
    assert sup_d[0] <= sup_e[0] <= chi0[0] <= chi0[1] <= sup_e[1] <= sup_d[1]
    save(fig, out, "pinning-the-ceiling.svg")
    return info


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        sys.exit(__doc__)
    root = args[0]
    here = os.path.dirname(os.path.abspath(__file__))
    out = here if os.path.basename(here) == "figures" else os.path.join(here, "figures")
    os.makedirs(out, exist_ok=True)
    info = dict(ceiling=fig_ceiling(root, out), windows=fig_windows(root, out))
    if "--info" in sys.argv:
        print(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
