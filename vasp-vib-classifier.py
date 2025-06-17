#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
from ase.io import read
from sklearn.decomposition import PCA
from collections import defaultdict


def load_vibmodes_from_outcar(inf='OUTCAR', exclude_imag=False):
    out = [line for line in open(inf) if line.strip()]
    ln = len(out)
    for line in out:
        if "NIONS =" in line:
            nions = int(line.split()[-1])
            break

    THz_index = []
    for ii in range(ln - 1, 0, -1):
        if '2PiTHz' in out[ii]:
            THz_index.append(ii)
        if 'Eigenvectors and eigenvalues of the dynamical matrix' in out[ii]:
            i_index = ii + 2
            break
    j_index = THz_index[0] + nions + 2

    freq_lines = [line for line in out[i_index:j_index] if '2PiTHz' in line]
    real_freq = [('f/i=' not in line) for line in freq_lines]
    omegas = [float(line.split()[3]) for line in freq_lines]  # 3rd is the real THz

    modes = [line.split()[3:6] for line in out[i_index:j_index]
             if ('dx' not in line) and ('2PiTHz' not in line)]
    omegas = np.array(omegas, dtype=float)
    modes = np.array(modes, dtype=float).reshape((-1, nions, 3))

    if exclude_imag:
        omegas = omegas[real_freq]
        modes = modes[real_freq]
        real_freq = [True for _ in omegas]

    return omegas, modes, real_freq


def classify_mode_species(mode, species):
    contributions = defaultdict(float)
    for atom_idx, vec in enumerate(mode):
        amp = np.linalg.norm(vec)
        contributions[species[atom_idx]] += amp ** 2

    total = sum(contributions.values())
    norm_contrib = {k: v / total for k, v in contributions.items()}

    sorted_species = sorted(norm_contrib.items(), key=lambda x: -x[1])
    if sorted_species[0][1] > 0.7:
        return f"{sorted_species[0][0]}-dominated"
    else:
        return "Mixed"


def classify_mode_localization(mode):
    # Flatten displacements and apply PCA
    X = mode.reshape(-1, 3)
    pca = PCA(n_components=1)
    pca.fit(X)
    explained = pca.explained_variance_ratio_[0]
    return "Collective" if explained > 0.8 else "Localized"


def write_xsf(imode, atoms, vector, scale=1.0,
              frequency=None, is_real=True,
              mode_type="", localization=""):
    vector = np.asarray(vector, dtype=float) * scale
    assert vector.shape == atoms.positions.shape
    pos_vec = np.hstack((atoms.positions, vector))
    nions = pos_vec.shape[0]
    chem_symbs = atoms.get_chemical_symbols()

    with open(f'mode_{imode:04d}.xsf', 'w') as out:
        header = ""
        if frequency is not None:
            freq_str = f"{abs(frequency):.6f}"
            tag = "Real" if is_real else "Imaginary"
            header += f"# Mode {imode}: {freq_str} cm-1 ({tag}, {mode_type}, {localization})\n"
        header += "CRYSTAL\nPRIMVEC\n"
        header += '\n'.join([
            ' '.join([f'{a:21.16f}' for a in vec])
            for vec in atoms.cell
        ])
        header += f"\nPRIMCOORD\n{nions:3d} 1\n"
        header += '\n'.join([
            f'{chem_symbs[i]:3s}' + ' '.join([f'{x:21.16f}' for x in pos_vec[i]])
            for i in range(nions)
        ])
        out.write(header)


def parse_args(cml):
    arg = argparse.ArgumentParser()
    arg.add_argument('-i', '--outcar', default='OUTCAR')
    arg.add_argument('-p', '--poscar', default='POSCAR')
    arg.add_argument('-m', '--mode', type=int, default=0,
                     help="0 = all modes, or specify mode index (1-based)")
    arg.add_argument('-s', '--scale', type=float, default=1.0)
    return arg.parse_args(cml)


def main(cml):
    p = parse_args(cml)
    atoms = read(p.poscar, format='vasp')
    species = atoms.get_chemical_symbols()

    if not all(os.path.exists(f) for f in ['OMEGAS.npy', 'MODES.npy', 'IS_REAL.npy']):
        omegas, modes, real_flags = load_vibmodes_from_outcar(p.outcar)
        modes /= np.sqrt(atoms.get_masses()[None, :, None])
        np.save('OMEGAS', omegas)
        np.save('MODES', modes)
        np.save('IS_REAL', real_flags)
    else:
        omegas = np.load('OMEGAS.npy')
        modes = np.load('MODES.npy')
        real_flags = np.load('IS_REAL.npy')

    n_modes = len(omegas)
    assert 0 <= p.mode <= n_modes

    if p.mode == 0:
        for i in range(n_modes):
            mode_type = classify_mode_species(modes[i], species)
            loc_type = classify_mode_localization(modes[i])
            write_xsf(i + 1, atoms, modes[i], p.scale,
                      frequency=omegas[i],
                      is_real=real_flags[i],
                      mode_type=mode_type,
                      localization=loc_type)
    else:
        i = p.mode - 1
        mode_type = classify_mode_species(modes[i], species)
        loc_type = classify_mode_localization(modes[i])
        write_xsf(p.mode, atoms, modes[i], p.scale,
                  frequency=omegas[i],
                  is_real=real_flags[i],
                  mode_type=mode_type,
                  localization=loc_type)


if __name__ == "__main__":
    main(sys.argv[1:])

