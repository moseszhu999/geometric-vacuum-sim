# Geometric Vacuum Simulation

**Simulating the Vacuum as a Born-Infeld Superfluid**

This repository contains the Python simulation code supporting the research on **Geometric Dynamics** and the **Vacuum Compressibility** hypothesis. These scripts were used to generate the data and plots for the paper: *"LHC Anomalies as Hydrodynamic Shockwaves in the Cosmic Superfluid"* (Zhu, 2025).

## Overview

The Institute of Geometric Dynamics (IGD) proposes that the vacuum is not empty space but a **Superfluid Bose-Einstein Condensate** composed of discrete "Space-Time Pixels" (Planck scale). This framework offers geometric explanations for several fundamental anomalies:

1.  **Fine Structure Constant ($\alpha \approx 1/137$)**: Arises from the geometric impedance mismatch between a toroidal vortex (electron) and the discrete vacuum grid.
2.  **LHC W-Boson Anomaly**: The $7\sigma$ mass shift observed by CDF II is modeled as a hydrodynamic shockwave in the compressible vacuum.

## Files

### 1. `vacuum_sim.py`
**The "Space-Time Pixel" Simulation**
- Models the vacuum as a 3D discrete lattice ($100^3$ voxels).
- Simulates an electron as a toroidal vortex ring.
- Calculates:
    - **Geometric Impedance ($Z_{eff}$)**: The energy cost to sustain a unit of topological flux.
    - **Reflection Coefficient ($R$)**: The probability of wave reflection from the vortex core (related to coupling efficiency $\alpha$).

### 2. `lhc_anomaly_plot.py`
**The Vacuum Shockwave Model**
- Models the W-boson mass shift as a function of collision energy $\sqrt{s}$.
- Uses a Rankine-Hugoniot saturation curve with a critical scale $\Lambda_{crit} \approx 246$ GeV (Higgs VEV).
- **Key Result**: Predicts the ~76 MeV mass shift at the LHC (13.6 TeV) while remaining consistent with LEP (209 GeV) data.
- Generates the prediction plot `lhc_anomaly_prediction.png`.

## Usage

### Prerequisites
- Python 3.x
- NumPy
- Matplotlib

### Installation
```bash
git clone https://github.com/moseszhu999/geometric-vacuum-sim.git
cd geometric-vacuum-sim
pip install -r requirements.txt
```

### Running the Simulations
To calculate vacuum impedance:
```bash
python vacuum_sim.py
```

To generate the LHC anomaly prediction plot:
```bash
python lhc_anomaly_plot.py
```

## Citation
If you use this code or the theoretical framework, please cite:

> Zhu, D. (2025). LHC Anomalies as Hydrodynamic Shockwaves in the Cosmic Superfluid. Zenodo. https://doi.org/10.5281/zenodo.17695129

> Zhu, D. (2025). Geometric Origin of the Fine Structure Constant via Vacuum Pixelation Efficiency. Zenodo. https://doi.org/10.5281/zenodo.17695341

## License
MIT License
