# VASP Vibrational Mode Classifier

A Python script to extract, classify, and export vibrational modes from VASP's `OUTCAR` and `POSCAR` files.

This tool:
- Parses vibrational modes from `OUTCAR`
- Classifies modes by **atomic species contribution** and **localization** (via PCA)
- Writes each mode as an `.xsf` file for visualization (e.g. in VESTA or VMD)


## Features

- Detects real vs. imaginary vibrational modes
- Species-dominance classification (e.g., "O-dominated", "Mixed")
- Localization analysis using principal component analysis
- Outputs `.xsf` files with annotated metadata


## Usage

```bash
python vasp-vib-classifier.py -i OUTCAR -p POSCAR
```

### Optional arguments:

| Flag | Description                              |
| ---- | ---------------------------------------- |
| `-i` | Path to `OUTCAR` (default: `OUTCAR`)     |
| `-p` | Path to `POSCAR` (default: `POSCAR`)     |
| `-m` | Mode index to export (default: all)      |
| `-s` | Displacement scale factor (default: 1.0) |

Example:

```bash
python vasp-vib-classifier.py -m 5 -s 0.5
```


## Output

The script generates `.xsf` files:

```
mode_0001.xsf, mode_0002.xsf, ...
```

Each file includes:

* Atomic coordinates
* Vibrational displacements (scaled)
* Frequency metadata
* Classification tags: species and localization


## Acknowledgments

This script is based on and extends functionality from:

**[QijingZheng/VaspVib2XSF](https://github.com/QijingZheng/VaspVib2XSF/blob/master/vasp2xsf.py)**
A Python tool for exporting VASP vibrational modes to XSF format.
