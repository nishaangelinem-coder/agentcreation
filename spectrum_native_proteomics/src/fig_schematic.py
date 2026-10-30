"""Fig. 1 (paradigm) and Fig. 3 (architectures) schematics."""
from figstyle import *
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


def box(ax, x, y, w, h, text, fc, ec=None, fs=6.3, tc=INK, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.012,rounding_size=0.015',
                                fc=fc, ec=ec or fc, lw=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs, color=tc,
            fontweight='bold' if bold else 'normal', wrap=True)


def arrow(ax, x0, y0, x1, y1, c=INK2, ls='-'):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=9, color=c, lw=1.1, ls=ls))


def fig1(path):
    fig, ax = plt.subplots(figsize=(7.2, 3.6)); ax.set_xlim(-0.02, 1.02); ax.set_ylim(0, 1); ax.axis('off')
    ax.text(0.01, 0.95, 'a  Identification-first paradigm', fontsize=9, fontweight='bold')
    ax.text(0.01, 0.46, 'b  Spectrum-native paradigm', fontsize=9, fontweight='bold')
    grey, lb = '#eeedea', '#dbe8f8'
    steps = ['LC-MS/MS\nraw MS1 + MS2', 'Database search\n(target-decoy, 1% FDR)', 'Peptide IDs\n(library-restricted)',
             'Protein quant.\n+ imputation', 'Statistics /\nclassifier']
    for i, s in enumerate(steps):
        box(ax, 0.01 + i * 0.198, 0.68, 0.165, 0.17, s, grey)
        if i: arrow(ax, 0.01 + i * 0.198 - 0.031, 0.765, 0.01 + i * 0.198 - 0.002, 0.765)
    ax.text(0.30, 0.6, 'loss: unidentified spectra, variants, PTMs, non-canonical peptides',
            fontsize=7, color=C[7], ha='left', style='italic')
    ax.text(0.30, 0.555, 'loss: stochastic DDA identification → missing values', fontsize=7, color=C[7], style='italic')
    box(ax, 0.01, 0.12, 0.165, 0.22, 'LC-MS/MS\nraw MS1 + MS2', lb)
    box(ax, 0.23, 0.24, 0.21, 0.13, 'Self-supervised pretraining\non unlabeled corpora\n(MassNet-style)', '#e5f4ec')
    box(ax, 0.23, 0.05, 0.21, 0.15, 'Spectrum encoder\n(transformer over peaks\nor MS1 feature tokens)', lb)
    arrow(ax, 0.335, 0.24, 0.335, 0.205, c=C[2])
    arrow(ax, 0.178, 0.2, 0.226, 0.13)
    outs = [('De novo sequencing\n(no database)', C[0]), ('Spectrum-to-phenotype\n(no identification)', C[1]),
            ('Representations / retrieval\n("proteomics GPT")', C[6])]
    for k, (s, c) in enumerate(outs):
        y = 0.30 - k * 0.125
        box(ax, 0.55, y - 0.005, 0.2, 0.1, s, 'white', ec=c, fs=6.0)
        arrow(ax, 0.443, 0.125, 0.547, y + 0.045, c=c)
    box(ax, 0.80, 0.07, 0.19, 0.26, 'Outputs\n• novel / variant peptides\n• cohort classification\n• attribution maps\n  at feature level', '#f6f5f2', fs=6.0)
    for k in range(3): arrow(ax, 0.752, 0.345 - k * 0.125, 0.797, 0.2, c=INK2)
    fig.savefig(path); plt.close(fig)


def fig3(path):
    fig, ax = plt.subplots(figsize=(7.2, 3.3)); ax.set_xlim(-0.02, 1.02); ax.set_ylim(0, 1); ax.axis('off')
    ax.text(0.0, 0.96, 'a  SpecNovo (MS2 → peptide)', fontsize=9, fontweight='bold')
    ax.text(0.53, 0.96, 'b  FeatureFormer (MS1 → phenotype)', fontsize=9, fontweight='bold')
    lb, gr, pk, yl = '#dbe8f8', '#e5f4ec', '#fbe3d7', '#fdf0cf'
    box(ax, 0.0, 0.70, 0.22, 0.16, 'Peaks {(m/z, I)} ≤ 80\n+ complementary m/z\n(M+2H−m/z)', '#f6f5f2')
    box(ax, 0.0, 0.47, 0.22, 0.15, 'Sinusoidal mass encoding\nλ ∈ [10⁻³, 10⁴] Da\n+ intensity projection', lb)
    box(ax, 0.0, 0.25, 0.22, 0.15, 'Transformer encoder\n3 layers, d=128, 4 heads', lb)
    box(ax, 0.26, 0.47, 0.22, 0.15, 'Decoder tokens: residue\n+ prefix mass + residual\nmass (precursor M, z)', pk)
    box(ax, 0.26, 0.25, 0.22, 0.15, 'Transformer decoder\n3 layers, cross-attention', pk)
    box(ax, 0.26, 0.03, 0.22, 0.15, 'Ion-evidence term\nI(b, y, b−H₂O, y−H₂O,\nb²⁺, y²⁺) → MLP → logit', yl)
    box(ax, 0.0, 0.03, 0.22, 0.15, 'Mass-constrained beam\nsearch (K=5, |ΔM|<0.05 Da)', '#f6f5f2')
    for (x0, y0, x1, y1) in [(0.11, 0.70, 0.11, 0.625), (0.11, 0.47, 0.11, 0.405), (0.22, 0.325, 0.257, 0.325),
                             (0.37, 0.47, 0.37, 0.405), (0.37, 0.18, 0.37, 0.247), (0.26, 0.1, 0.223, 0.1),
                             (0.26, 0.29, 0.2, 0.18)]:
        arrow(ax, x0, y0, x1, y1)
    box(ax, 0.53, 0.70, 0.21, 0.16, 'MS1 feature detection\n(m/z, RT, log₂ I)\nno MS2 identification', '#f6f5f2')
    box(ax, 0.53, 0.47, 0.21, 0.15, 'ID-free alignment to\nconsensus map (RT-shift\n+ 15 ppm / 0.7 min match)', lb)
    box(ax, 0.53, 0.25, 0.21, 0.15, 'Token = feature emb.\n+ sin(m/z), sin(RT)\n+ value | missing token', lb)
    box(ax, 0.78, 0.47, 0.21, 0.15, 'SSL: masked intensity\nmodelling (30 %)\non unlabeled runs', gr)
    box(ax, 0.78, 0.25, 0.21, 0.15, 'Transformer encoder\n2 layers, d=64, [CLS]', lb)
    box(ax, 0.78, 0.03, 0.21, 0.15, 'Head: CLS logit +\nadditive token evidence\n→ P(case)', pk)
    box(ax, 0.53, 0.03, 0.21, 0.15, 'Gradient × input\nattribution per feature\n(biomarker discovery)', yl)
    for (x0, y0, x1, y1) in [(0.635, 0.70, 0.635, 0.625), (0.635, 0.47, 0.635, 0.405), (0.74, 0.325, 0.777, 0.325),
                             (0.885, 0.47, 0.885, 0.405), (0.885, 0.25, 0.885, 0.185), (0.78, 0.1, 0.743, 0.1)]:
        arrow(ax, x0, y0, x1, y1)
    fig.savefig(path); plt.close(fig)


if __name__ == '__main__':
    fig1('figures/fig1_paradigm.png'); fig3('figures/fig3_architecture.png')
