#!/usr/bin/env bash
# 编译：生成图 2（python3 + matplotlib）-> xelatex -> bibtex -> xelatex -> xelatex
set -e
cd "$(dirname "$0")"
mkdir -p figures
python3 - <<'PY'
# 图 2：命题 1 在基线参数下的数值示意，全部由式 (4)、式 (A.4) 计算
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

for name in ['Noto Serif CJK SC', 'Source Han Serif SC', 'SimSun', 'Songti SC']:
    if any(f.name == name for f in font_manager.fontManager.ttflist):
        plt.rcParams['font.family'] = [name]
        break
plt.rcParams.update({'mathtext.fontset': 'stix', 'font.size': 9, 'axes.unicode_minus': False,
                     'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6})

Rbar, vR, k, rho = 1.15, 0.04, 0.35, 0.25
U = Rbar ** 2 + vR
A = Rbar + rho * U / k
mx = A / vR
m = np.logspace(np.log10(5), np.log10(1000), 400)

fig, (a, b) = plt.subplots(1, 2, figsize=(6.4, 2.55))
a.loglog(m, vR + A / m, color='#1f3b73', lw=1.6, label=r'$V(m)=v_R+A/m$')
a.loglog(m, A / m, color='#c0504d', lw=1.0, ls='--', label=r'抽样项 $A/m$')
a.loglog(m, np.full_like(m, vR), color='#4f8a3c', lw=1.0, ls=':', label=r'环境项 $v_R$')
a.axvline(mx, color='0.45', lw=0.7)
a.annotate(r'$m_{\times}=%.2f$' % mx, (mx, 0.3), xytext=(6, 0), textcoords='offset points',
           fontsize=8.5, color='0.25')
a.set_xlabel(r'病例量 $m$')
a.set_ylabel('误差方差')
a.set_title('(a) 误差方差及其两部分', fontsize=9)
a.legend(frameon=False, fontsize=8, loc='lower left')

g = 1 / (2 * (1 + m / mx))
b.semilogx(m, g, color='#1f3b73', lw=1.6)
for mm, lab in ((10, r'$m=10$'), (mx, r'$m=m_{\times}$'), (200, r'$m=200$')):
    y = 1 / (2 * (1 + mm / mx))
    b.plot([mm], [y], 'o', ms=3.5, color='#c0504d')
    b.annotate('%s: %.3f' % (lab, y), (mm, y), xytext=(5, 4), textcoords='offset points', fontsize=8)
b.set_ylim(0, 0.5)
b.set_xlabel(r'病例量 $m$')
b.set_ylabel(r'加倍 $m$ 的相对降幅')
b.set_title(r'(b) 式 (A.4)：$1/[2(1+m/m_{\times})]$', fontsize=9)
for ax in (a, b):
    ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout(w_pad=2.0)
fig.savefig('figures/fig-prop1.png', dpi=400)
PY
xelatex -interaction=nonstopmode -halt-on-error theory-derivation
bibtex theory-derivation
xelatex -interaction=nonstopmode -halt-on-error theory-derivation
xelatex -interaction=nonstopmode -halt-on-error theory-derivation
