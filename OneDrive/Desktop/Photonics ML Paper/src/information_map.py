"""
Task 6: information map.

With only 15 raw structures living natively in a 3D parameter space,
t-SNE/UMAP are not appropriate here -- both are designed to unfold LOCAL
neighborhood structure of high-dimensional, large-n data, and are known to
produce unstable, easily over-interpreted embeddings for n this small in a
space this low-dimensional. PCA on 3D data is just "find the best-fit
plane," which is exactly what's needed to visualize occupancy without
manufacturing structure that isn't there -- so PCA is used here, and
UMAP/t-SNE are deliberately skipped (a documented methodological choice,
not an oversight).

Answers "do current process monitors occupy only a small region of
Fisher-information space" using the SAME numeric evidence as atlas.py /
atlas_vs_fisher.py (within/cross-class angle stats, empty-wedge analysis)
-- the PCA plot illustrates the conclusion, it doesn't establish it.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from atlas import load_records, label


def pca_2d(X):
    """From-scratch PCA via SVD (no sklearn dependency)."""
    mean = X.mean(axis=0)
    Xc = X - mean
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    explained = (S**2) / np.sum(S**2)
    proj = Xc @ Vt.T
    return proj, explained, Vt


if __name__ == "__main__":
    records = load_records()
    labels = [label(r) for r in records]
    scaled = np.array([r['sens_vector_scaled'] for r in records])
    families = [f"{r['geometry']}-{r['polarization']}" for r in records]

    proj, explained, components = pca_2d(scaled)
    print("=== PCA of fab-tolerance-scaled sensitivity vectors (3D -> 2D) ===")
    print(f"Variance explained: PC1={explained[0]*100:.1f}%  PC2={explained[1]*100:.1f}%  "
          f"PC3={explained[2]*100:.1f}%  (PC1+PC2={100*(explained[0]+explained[1]):.1f}%)")
    print(f"PC1 loadings (width, thickness, angle): {components[0].round(3)}")
    print(f"PC2 loadings (width, thickness, angle): {components[1].round(3)}")

    print("\nPer-family centroid in PC space:")
    for fam in sorted(set(families)):
        idx = [i for i, f in enumerate(families) if f == fam]
        c = proj[idx].mean(axis=0)
        print(f"  {fam:>10}: PC1={c[0]:+.4f}  PC2={c[1]:+.4f}  (n={len(idx)})")

    # inter-centroid separation vs intra-family spread, quantified
    fam_set = sorted(set(families))
    centroids = {fam: proj[[i for i, f in enumerate(families) if f == fam]].mean(axis=0) for fam in fam_set}
    spreads = {fam: np.std(proj[[i for i, f in enumerate(families) if f == fam]], axis=0).mean() for fam in fam_set}
    print("\nInter-family centroid distances vs intra-family spread (occupancy check):")
    for i, f1 in enumerate(fam_set):
        for f2 in fam_set[i+1:]:
            d = np.linalg.norm(centroids[f1] - centroids[f2])
            avg_spread = (spreads[f1] + spreads[f2]) / 2
            print(f"  {f1} <-> {f2}: centroid dist={d:.4f}, avg intra-family spread={avg_spread:.4f}, "
                  f"ratio={d/avg_spread if avg_spread>0 else float('inf'):.1f}x")

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = {'strip-TE': 'tab:blue', 'strip-TM': 'tab:orange', 'slot-TE': 'tab:green',
              'slot-TM': 'tab:olive', 'rib-TE': 'tab:red', 'swg-TE': 'tab:purple'}
    for fam in fam_set:
        idx = [i for i, f in enumerate(families) if f == fam]
        ax.scatter(proj[idx, 0], proj[idx, 1], label=fam, color=colors.get(fam, 'gray'), s=60)
        for i in idx:
            ax.annotate(labels[i].split('-')[-1], (proj[i, 0], proj[i, 1]), fontsize=7,
                        xytext=(4, 4), textcoords='offset points')
    ax.set_xlabel(f'PC1 ({explained[0]*100:.0f}% var)')
    ax.set_ylabel(f'PC2 ({explained[1]*100:.0f}% var)')
    ax.set_title('Information map: PCA of fab-tolerance-scaled sensitivity vectors\n'
                 f'(n={len(records)} structures)', fontsize=11)
    ax.legend()
    ax.axhline(0, color='0.85', lw=0.7)
    ax.axvline(0, color='0.85', lw=0.7)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig('figures/fig2_information_map.png', dpi=150)
    plt.close(fig)
    print("\nFigure saved: figures/fig2_information_map.png")
