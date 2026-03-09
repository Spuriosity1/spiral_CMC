#!/usr/bin/env python3
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def wrap_to_first_bz(k_points, B):
    inv_B = np.linalg.inv(B)
    frac_k = k_points @ inv_B
    frac_k_wrapped = frac_k - np.round(frac_k)

    n0 = k_points.shape[0]
    n1 = k_points.shape[1]
    n2k = k_points.shape[2]
    n2 = (n2k-1)*2

    return frac_k_wrapped @ B

def load_data(ssf_path, bz_path):
    data = {'observables': {}}
#    bz_path = list(Path(ssf_path).parent.glob("*.BZ.h5"))[0]
    
    with h5py.File(ssf_path, 'r') as f:
        data['n_samples'] = f['n_samples'][0] if 'n_samples' in f else 1.0
        for obs in f['ssf'].keys():
            data['observables'][obs] = f['ssf'][obs][:]

    with h5py.File(bz_path, 'r') as f:
        k_pts = f['k_points'][:]
        inv_prim = f['inverse_primitive'][:]

        B = 2 * np.pi * inv_prim
        data['k_points'] = wrap_to_first_bz(k_pts, B)
        
    return data, B

# -------------------------------------------------------------------
# Select points whose projection onto the normal vector is ≈ projection
# -------------------------------------------------------------------
def take_cut(k_data, normal_vector, projection, tol=1e-4):
    normal_vector = normal_vector / np.linalg.norm(normal_vector)
    proj = k_data @ normal_vector
    mask = np.abs(proj - projection) < tol
    return mask, proj


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Path to .ssf file')
    parser.add_argument('bz_file', help='Path to .BZ file')
    parser.add_argument('--observable', help='e.g., Szz', default='SdotS')
    parser.add_argument('--plane', nargs=3, type=float, default=[1,0,0],
                        help='normal vector for the plane to be plotted')
    parser.add_argument('--coord', type=float, default=0.0,
                        help='Projection value onto plane normal')
    parser.add_argument('--xlim', type=float, default=10)
    parser.add_argument('--ylim', type=float, default=10)
    parser.add_argument('--repeat', type=int, nargs=3, default=[0,0,0])
    args = parser.parse_args()

    data, B = load_data(args.file, args.bz_file)
    kpts = data['k_points']
    obs = data['observables'][args.observable]

    # ------------------------------------------------------------
    # 1. Take the cut
    # ------------------------------------------------------------
    mask, proj_vals = take_cut(kpts, np.array(args.plane), args.coord)
    if(obs.shape != mask.shape):
        print("Incorrect BZ file supplied: observable {obs.shape} != BZ {mask.shape}")

    k_cut = kpts[mask,:]
    s_cut = obs[mask]

    print(f"{sum(mask.ravel())} points in cut")

    print(k_cut)

    # ------------------------------------------------------------
    # 2. Build a 2D coordinate system inside the plane
    # ------------------------------------------------------------
    n = np.array(args.plane, dtype=float)
    n /= np.linalg.norm(n)

    # pick a vector not parallel to n
    tmp = np.array([1,0,0]) if abs(n[0]) < 0.9 else np.array([0,1,0])
    e1 = np.cross(n, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)

    Emat = np.array([e1, e2]).T

    # project k-points into plane basis
    c1 = k_cut @ e1
    c2 = k_cut @ e2

    # project basis vectors into plane basis
    Bproj = B @ Emat

    # ------------------------------------------------------------
    # 4. Plot
    # ------------------------------------------------------------
    plt.figure(figsize=(9, 7))

    repeat = args.repeat

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')

    xlim = args.xlim
    ylim = args.ylim
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-ylim, ylim)
    plt.xlabel("k₁ (Å⁻¹)")
    plt.ylabel("k₂ (Å⁻¹)")
    plt.title(f"{args.observable} cut normal to {args.plane}")
    plt.tight_layout()

    # Estimate marker size from typical point spacing
    diffs = np.diff(np.sort(c1), prepend=np.min(c1)-1)
    spacing = np.median(diffs[diffs > 1e-6])
    # Convert data spacing to points: figure size * dpi / data range
    px_per_unit = min(9 * 72 / (2 * xlim), 7 * 72 / (2 * ylim))
    marker_size = (spacing * px_per_unit) ** 2

    for i0 in range(-repeat[0], repeat[0]+1):
        for i1 in range(-repeat[1], repeat[1]+1):
            for i2 in range(-repeat[2], repeat[2]+1):
                c1_q, c2_q = np.transpose(
                        np.transpose([c1, c2]) + [i0,i1,i2]@Bproj )
                mesh = ax.scatter( c1_q,  c2_q, c=s_cut / data['n_samples'],
                                   cmap='viridis', s=marker_size,
                                   marker='s', edgecolors=None)
                ax.scatter(-c1_q, -c2_q, c=s_cut / data['n_samples'],
                                   cmap='viridis', s=marker_size,
                                   marker='s', edgecolors=None)

    plt.colorbar(mesh, label='S(q)')
    plt.show()

if __name__ == "__main__":
    main()

