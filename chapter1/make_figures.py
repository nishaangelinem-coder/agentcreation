import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import sys
out = sys.argv[1]
plt.rcParams["font.family"] = "DejaVu Sans"

def box(ax, x, y, w, h, text, fc, ec="#333333", fs=9, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=fc, ec=ec, lw=1.1))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", wrap=True)

def arrow(ax, x1, y1, x2, y2, color="#444444"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.3, color=color))

# Figure 1: evolution of agriculture 1.0 -> 5.0
fig, ax = plt.subplots(figsize=(10, 3.6)); ax.set_xlim(0, 10); ax.set_ylim(0, 3.6); ax.axis("off")
stages = [("Agriculture 1.0", "Manual labour and\nanimal power", "#e8f5e9"),
          ("Agriculture 2.0", "Mechanisation,\nsynthetic inputs\n(Green Revolution)", "#c8e6c9"),
          ("Agriculture 3.0", "GNSS guidance,\nvariable-rate tech,\nprecision farming", "#a5d6a7"),
          ("Agriculture 4.0", "IoT, big data,\ncloud, machine\nlearning", "#81c784"),
          ("Agriculture 5.0", "AI-driven autonomy,\nrobotics, digital\ntwins, human-centric", "#66bb6a")]
w = 1.7; gap = 0.25
for i, (t, d, c) in enumerate(stages):
    x = 0.15 + i * (w + gap)
    box(ax, x, 2.45, w, 0.75, t, c, fs=10, bold=True)
    box(ax, x, 0.75, w, 1.45, d, "#ffffff", fs=8.5)
    if i < len(stages) - 1:
        arrow(ax, x + w + 0.02, 2.82, x + w + gap - 0.02, 2.82)
arrow(ax, 0.15, 0.35, 9.8, 0.35, "#2e7d32")
ax.text(5, 0.08, "Increasing data intensity, connectivity and decision autonomy", ha="center", fontsize=9, style="italic")
fig.savefig(f"{out}/Fig1_evolution.png", dpi=300, bbox_inches="tight"); plt.close(fig)

# Figure 2: layered architecture
fig, ax = plt.subplots(figsize=(8, 6.2)); ax.set_xlim(0, 8); ax.set_ylim(0, 6.6); ax.axis("off")
layers = [("Layer 1: Sensing and Data Acquisition", "Soil & weather sensors, cameras, multispectral/hyperspectral imagers,\nUAVs, satellites, machine telemetry, farmer records", "#e3f2fd"),
          ("Layer 2: Connectivity and Edge", "LoRaWAN, NB-IoT, Wi-Fi, 4G/5G, edge gateways, embedded/TinyML devices", "#bbdefb"),
          ("Layer 3: Data Management", "Cleaning, fusion, geo-referencing, annotation, cloud storage, data governance", "#90caf9"),
          ("Layer 4: AI and Analytics", "ML/DL models, computer vision, forecasting, optimisation, knowledge graphs", "#64b5f6"),
          ("Layer 5: Decision Support and Actuation", "Advisories, prescription maps, variable-rate machinery, robots, irrigation control", "#42a5f5")]
for i, (t, d, c) in enumerate(layers):
    y = 0.3 + i * 1.25
    box(ax, 0.3, y, 6.4, 1.0, "", c)
    ax.text(3.5, y + 0.7, t, ha="center", fontsize=10, fontweight="bold")
    ax.text(3.5, y + 0.3, d, ha="center", va="center", fontsize=8)
    if i < len(layers) - 1:
        arrow(ax, 3.5, y + 1.0, 3.5, y + 1.25)
box(ax, 6.95, 0.3, 0.8, 6.0, "", "#fff3e0")
ax.text(7.35, 3.3, "Security, privacy, explainability & ethics", rotation=90, ha="center", va="center", fontsize=9, fontweight="bold")
fig.savefig(f"{out}/Fig2_architecture.png", dpi=300, bbox_inches="tight"); plt.close(fig)

# Figure 3: taxonomy of AI techniques mapped to tasks
fig, ax = plt.subplots(figsize=(10, 5.4)); ax.set_xlim(0, 10); ax.set_ylim(0, 5.6); ax.axis("off")
box(ax, 3.6, 4.65, 2.8, 0.7, "Artificial Intelligence\nin Agriculture", "#fff9c4", fs=10, bold=True)
cats = [("Knowledge-based\nsystems", "Expert systems,\nfuzzy logic,\nontologies", "Advisory services,\npest diagnosis rules"),
        ("Classical machine\nlearning", "SVM, random forest,\ngradient boosting,\nk-means", "Yield prediction,\nsoil property mapping"),
        ("Deep learning", "CNN, RNN/LSTM,\nTransformers", "Disease detection,\nweed/fruit detection"),
        ("Reinforcement\nlearning", "Q-learning,\npolicy gradients", "Irrigation scheduling,\nrobot navigation"),
        ("Foundation &\ngenerative models", "Vision/language\nfoundation models", "Farmer chatbots,\nsynthetic data")]
w = 1.78; gap = 0.2
for i, (c, m, a) in enumerate(cats):
    x = 0.1 + i * (w + gap)
    box(ax, x, 3.2, w, 0.85, c, "#e1bee7", fs=9, bold=True)
    box(ax, x, 1.85, w, 1.05, m, "#f3e5f5", fs=8.3)
    box(ax, x, 0.35, w, 1.05, a, "#e8f5e9", fs=8.3)
    arrow(ax, 5.0, 4.65, x + w/2, 4.05)
    arrow(ax, x + w/2, 3.2, x + w/2, 2.9)
    arrow(ax, x + w/2, 1.85, x + w/2, 1.4)
ax.text(10.05, 2.37, "Methods", fontsize=8.5, style="italic", color="#6a1b9a", rotation=90, va="center")
ax.text(10.05, 0.87, "Tasks", fontsize=8.5, style="italic", color="#1b5e20", rotation=90, va="center")
fig.savefig(f"{out}/Fig3_taxonomy.png", dpi=300, bbox_inches="tight"); plt.close(fig)

# Figure 4: closed-loop data-to-decision cycle
import numpy as np
fig, ax = plt.subplots(figsize=(7, 7)); ax.set_xlim(-3.6, 3.6); ax.set_ylim(-3.6, 3.6); ax.axis("off"); ax.set_aspect("equal")
steps = ["1. Observe\n(sense the field)", "2. Collect &\ntransmit data", "3. Pre-process\n& fuse data",
         "4. Learn & infer\n(AI models)", "5. Decide\n(recommendation)", "6. Act\n(machinery, robots,\nfarmer)", "7. Evaluate\n& feed back"]
n = len(steps); r = 2.55
cols = ["#c8e6c9", "#b2dfdb", "#b3e5fc", "#d1c4e9", "#ffe0b2", "#ffccbc", "#f0f4c3"]
pts = []
for i, s in enumerate(steps):
    th = np.pi/2 - 2*np.pi*i/n
    x, y = r*np.cos(th), r*np.sin(th); pts.append((x, y))
    box(ax, x-0.85, y-0.45, 1.7, 0.9, s, cols[i], fs=8.3)
for i in range(n):
    (x1, y1), (x2, y2) = pts[i], pts[(i+1) % n]
    arrow(ax, x1 + 0.36*(x2-x1), y1 + 0.36*(y2-y1), x1 + 0.64*(x2-x1), y1 + 0.64*(y2-y1), "#2e7d32")
ax.text(0, 0, "Continuous\nlearning loop", ha="center", va="center", fontsize=11, fontweight="bold", color="#1b5e20")
fig.savefig(f"{out}/Fig4_loop.png", dpi=300, bbox_inches="tight"); plt.close(fig)
print("done")
