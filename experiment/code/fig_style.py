"""Figure style shared by make_fig*.py: imitates Eales & Riley (2024) Epidemics 46:100742
(ggplot2 theme_bw look: white panel, light grey grid, thin black frame; ColorBrewer Dark2;
legends boxed inside the panel with a title; bold upper-case panel letters outside the top-left).

Fonts: matplotlib cannot resolve "Times New Roman" in this image, so Latin text uses Liberation
Serif (metric-compatible Times clone) and Chinese text uses Noto Serif CJK SC (Song/Ming class).
Every text artist is single-script; numbers in Chinese figures stay in the Latin font.
Collision checks reuse .claude/skills/experiment-suite/figure_examples/style_kit.py.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.text import Text

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / ".claude/skills/experiment-suite/figure_examples"))
from style_kit import (assert_no_legend_reference_line_collisions,  # noqa: E402
                       assert_no_text_collisions)

EXP = Path(__file__).resolve().parents[1]
FIGDIR = EXP / "figures"

LATIN = "Liberation Serif"
CJK_PATH = "/usr/share/fonts/opentype/noto/NotoSerifCJK-Regular.ttc"
font_manager.fontManager.addfont(CJK_PATH)
CJK = "Noto Serif CJK SC"
CJK_RE = re.compile(r"[　-〿㐀-鿿＀-￯]")
LAT_RE = re.compile(r"[A-Za-z0-9]")

DARK2 = ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"]
GREY = "#666666"

MM = 1 / 25.4
FULL_W = 174 * MM           # Elsevier double-column width
HALF_W = 84 * MM


def cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lang", choices=["en", "zh"], default="en")
    return ap.parse_args().lang


def setup(lang: str):
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": [LATIN],
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7,
        "legend.title_fontsize": 7,
        "axes.edgecolor": "black",
        "axes.linewidth": 0.6,
        "axes.grid": True,
        "grid.color": "#D9D9D9",
        "grid.linewidth": 0.4,
        "axes.facecolor": "white",
        "axes.axisbelow": True,
        "axes.prop_cycle": matplotlib.cycler(color=DARK2),
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.frameon": True,
        "legend.edgecolor": "#BBBBBB",
        "legend.fancybox": False,
        "legend.framealpha": 0.95,
        "lines.linewidth": 1.0,
        "mathtext.fontset": "custom",
        "mathtext.rm": LATIN,
        "mathtext.it": f"{LATIN}:italic",
        "mathtext.bf": f"{LATIN}:bold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
        "savefig.dpi": 300,
    })


def T(lang, en, zh):
    """Pick label by language."""
    return zh if lang == "zh" else en


def fp_for(text: str) -> FontProperties | None:
    return FontProperties(family=CJK) if CJK_RE.search(text) else None


def xlabel(ax, s):
    ax.set_xlabel(s, fontproperties=fp_for(s), fontsize=8)


def ylabel(ax, s):
    ax.set_ylabel(s, fontproperties=fp_for(s), fontsize=8)


def title(ax, s):
    ax.set_title(s, fontproperties=fp_for(s), fontsize=8, loc="left")


def text(ax, x, y, s, **kw):
    kw.setdefault("fontsize", 7)
    return ax.text(x, y, s, fontproperties=fp_for(s), **kw)


def legend(ax, title_s=None, **kw):
    """Boxed in-panel legend with a title; CJK entries/titles get the CJK font."""
    kw.setdefault("loc", "best")
    lg = ax.legend(title=title_s, **kw)
    for t in lg.get_texts():
        if CJK_RE.search(t.get_text()):
            t.set_fontproperties(FontProperties(family=CJK, size=7))
    if title_s and CJK_RE.search(title_s):
        lg.get_title().set_fontproperties(FontProperties(family=CJK, size=7))
    lg.get_frame().set_linewidth(0.5)
    return lg


def panel_letter(ax, letter, dx=-0.16, dy=1.03):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            ha="left", va="bottom", family="serif")


def check_fonts(fig, allow_bold={"A", "B", "C", "D", "E", "F"}):
    """Each visible text artist is single-script, >= 7 pt, resolves to Liberation Serif (Latin)
    or Noto Serif CJK SC (CJK) without fallback; bold only for panel letters."""
    fig.canvas.draw()
    for t in fig.findobj(match=Text):
        s = t.get_text()
        if not t.get_visible() or not s.strip():
            continue
        cjk, lat = bool(CJK_RE.search(s)), bool(LAT_RE.search(re.sub(r"\$[^$]*\$", "", s)))
        if cjk and lat:
            raise RuntimeError(f"mixed-script text artist: {s!r}")
        if not (cjk or lat or "$" in s):
            continue
        if t.get_fontsize() < 7:
            raise RuntimeError(f"font < 7 pt: {s!r}")
        path = font_manager.findfont(t.get_fontproperties(), fallback_to_default=False)
        name = FontProperties(fname=path).get_name()
        want = CJK if cjk else LATIN
        if name != want:
            raise RuntimeError(f"font {name!r} != {want!r} for {s!r}")
        w = t.get_fontweight()
        if (w == "bold" or (isinstance(w, int) and w > 500)) and s not in allow_bold:
            raise RuntimeError(f"bold text outside panel letters: {s!r}")


def save(fig, name: str, lang: str, checks=True):
    out = FIGDIR / lang
    out.mkdir(parents=True, exist_ok=True)
    if checks:
        check_fonts(fig)
        assert_no_text_collisions(fig)
        assert_no_legend_reference_line_collisions(fig)
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(out / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("saved", out / name)
