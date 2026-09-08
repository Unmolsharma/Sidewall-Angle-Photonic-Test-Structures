"""
Task 2: sensitivity atlas.

Treats every simulated structure's normalized (fab-tolerance-scaled)
sensitivity vector (dneff/dw, dneff/dt, dneff/da) as a point/direction in
3D fabrication-parameter space. Quantifies, for every pair of structures:
  - vector angle (degrees, between unit vectors)
  - cosine similarity
  - Euclidean distance (between the SCALED, not unit, vectors -- this is
    deliberately different information than angle/cosine: it also reflects
    magnitude differences, i.e. "how much information," not just direction)

Produces numeric matrices (CSV) plus a multi-panel figure (candidate
Figure 1): three 2D projections (w-t, w-a, t-a) as origin-anchored
direction arrows, and a hierarchically-clustered cosine-similarity heatmap.
No conclusions are read off the plot alone -- every visual claim below is
backed by a printed number from the matrices.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform

AXIS_NAMES = ['width', 'thickness', 'angle']
AXIS_PAIRS = [(0, 1, 'width vs thickness'), (0, 2, 'width vs sidewall angle'), (1, 2, 'thickness vs sidewall angle')]


def load_records(path='data/sensitivity_db.json'):
    with open(path) as f:
        return json.load(f)


def label(r):
    return f"{r['geometry']}-{r['polarization']}-{r['width_nm']:.0f}"


def build_matrices(records):
    labels = [label(r) for r in records]
    unit = np.array([r['sens_vector_unit'] for r in records])     # for angle/cosine
    scaled = np.array([r['sens_vector_scaled'] for r in records])  # for Euclidean distance
    n = len(records)

    cos_sim = unit @ unit.T
    cos_sim = np.clip(cos_sim, -1.0, 1.0)
    angle_deg = np.degrees(np.arccos(cos_sim))
    euclid = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            euclid[i, j] = np.linalg.norm(scaled[i] - scaled[j])

    return labels, unit, scaled, cos_sim, angle_deg, euclid


def save_matrix_csv(path, labels, M):
    with open(path, 'w') as f:
        f.write(',' + ','.join(labels) + '\n')
        for lab, row in zip(labels, M):
            f.write(lab + ',' + ','.join(f'{v:.6f}' for v in row) + '\n')


def numeric_summary(labels, records, angle_deg, cos_sim):
    """Quantify clustering/redundancy WITHOUT relying on the plot: for every
    pair, classify as same-class (same geometry+polarization) or
    cross-class, and report angle statistics for each group."""
    n = len(records)
    same_class, cross_class = [], []
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            same = (records[i]['geometry'] == records[j]['geometry'] and
                    records[i]['polarization'] == records[j]['polarization'])
            entry = (labels[i], labels[j], angle_deg[i, j], cos_sim[i, j])
            pairs.append(entry)
            (same_class if same else cross_class).append(entry)

    def stats(group, name):
        if not group:
            return f"{name}: (none)"
        angs = [g[2] for g in group]
        return (f"{name}: n={len(group)}  angle min/mean/max = "
                f"{min(angs):.1f}/{np.mean(angs):.1f}/{max(angs):.1f} deg")

    lines = [stats(same_class, "within-class pairs (same geometry+pol, different width)"),
             stats(cross_class, "cross-class pairs (different geometry or polarization)")]

    most_redundant = min(pairs, key=lambda p: p[2])
    most_orthogonal = max(pairs, key=lambda p: p[2])
    lines.append(f"most redundant pair (smallest angle): {most_redundant[0]} <-> {most_redundant[1]} "
                 f"= {most_redundant[2]:.2f} deg (cos={most_redundant[3]:.4f})")
    lines.append(f"most orthogonal pair (largest angle): {most_orthogonal[0]} <-> {most_orthogonal[1]} "
                 f"= {most_orthogonal[2]:.2f} deg (cos={most_orthogonal[3]:.4f})")
    return "\n".join(lines), pairs


def make_figure(labels, unit, cos_sim, records, out_path='figures/fig1_sensitivity_atlas.png'):
    fig = plt.figure(figsize=(14, 10))
    colors = {'strip-TE': 'tab:blue', 'strip-TM': 'tab:orange', 'slot-TE': 'tab:green',
              'slot-TM': 'tab:olive', 'rib-TE': 'tab:red', 'swg-TE': 'tab:purple'}

    def color_for(r):
        return colors.get(f"{r['geometry']}-{r['polarization']}", 'gray')

    # --- 3 projection panels ---
    for k, (ax_i, ax_j, title) in enumerate(AXIS_PAIRS):
        ax = fig.add_subplot(2, 2, k+1)
        for r, v in zip(records, unit):
            c = color_for(r)
            ax.annotate('', xy=(v[ax_i], v[ax_j]), xytext=(0, 0),
                        arrowprops=dict(arrowstyle='->', color=c, lw=1.6, alpha=0.85))
            ax.scatter([v[ax_i]], [v[ax_j]], s=8, color=c)
        ax.axhline(0, color='0.8', lw=0.7)
        ax.axvline(0, color='0.8', lw=0.7)
        ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1)
        ax.set_aspect('equal')
        ax.set_xlabel(AXIS_NAMES[ax_i]); ax.set_ylabel(AXIS_NAMES[ax_j])
        ax.set_title(title, fontsize=10)

    handles = [plt.Line2D([0], [0], color=c, lw=2, label=k) for k, c in colors.items()]
    fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.98, 0.98), fontsize=8)

    # --- clustered cosine-similarity heatmap ---
    ax = fig.add_subplot(2, 2, 4)
    dist = 1 - cos_sim
    np.fill_diagonal(dist, 0)
    condensed = squareform(np.clip(dist, 0, None), checks=False)
    Z = linkage(condensed, method='average')
    order = leaves_list(Z)
    reordered = cos_sim[np.ix_(order, order)]
    im = ax.imshow(reordered, vmin=-1, vmax=1, cmap='RdBu_r')
    ax.set_xticks(range(len(labels))); ax.set_xticklabels([labels[i] for i in order], rotation=90, fontsize=6)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels([labels[i] for i in order], fontsize=6)
    ax.set_title('cosine similarity (hierarchically clustered)', fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle('Sensitivity atlas: normalized fabrication-parameter direction vectors', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return order, Z


if __name__ == "__main__":
    import os
    os.makedirs('data', exist_ok=True)
    os.makedirs('figures', exist_ok=True)

    records = load_records()
    labels, unit, scaled, cos_sim, angle_deg, euclid = build_matrices(records)

    save_matrix_csv('data/atlas_cosine_similarity.csv', labels, cos_sim)
    save_matrix_csv('data/atlas_angle_deg.csv', labels, angle_deg)
    save_matrix_csv('data/atlas_euclidean_distance.csv', labels, euclid)

    summary, pairs = numeric_summary(labels, records, angle_deg, cos_sim)
    print("=== Sensitivity atlas: numeric summary ===")
    print(summary)

    order, Z = make_figure(labels, unit, cos_sim, records)
    print(f"\nFigure saved: figures/fig1_sensitivity_atlas.png")
    print(f"Clustering leaf order: {[labels[i] for i in order]}")

    with open('data/atlas_summary.txt', 'w') as f:
        f.write(summary + '\n')
