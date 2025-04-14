# 🌀 Spintronics in Terahertz Communication

This repository contains simulation code and numerical analysis tools related to spintronic device behavior in terahertz (THz) communication systems. The codebase is part of an MSc research project focusing on spin wave generation, signal amplification, and signal detection using spintronic materials and devices.

## 📚 Project Highlights

- Simulation of **magnetisation dynamics** in various materials (Al, Au, Ga, GaAs).
- Modelling of **spintronic devices**: Spin Valves, Magnetic Tunnel Junctions (MTJs), and Spin Torque Nano Oscillators (STNOs).
- Analysis of **spin-related effects**: spin pumping, spin wave generation, and spin transfer torque (STT).
- Comparison of **numerical integration methods**: RK4, AB2, and ode45 (RK45).

## 🧪 Features

- Implements the Landau-Lifshitz-Gilbert (LLG) equation with damping and torque terms.
- Time-domain simulation of dynamic magnetisation behavior in THz fields.
- Clean and modular MATLAB code structure with plotting support.
- Realistic modeling using experimental parameters.

## 📁 Repository Structure

```
Spintronics_Bristol/
├── simulation/
│   ├── llg_solver.m                # LLG equation solver
│   ├── simulate_spin_effects.m     # Spin pumping, spin wave, and STT simulation
│   ├── integration_comparison.m    # Comparison of RK4, AB2, and ode45
│   └── materials_params.m          # Magnetic parameters for materials
├── devices/
│   ├── simulate_devices.m          # STNO, MTJ, Spin Valve modeling
│   └── device_params.m             # Device-specific configurations
├── plots/                          # Auto-generated result figures
├── data/                           # Optional saved output files
└── README.md
```

## ⚙️ Requirements

- MATLAB R2022a or later
- No third-party toolboxes required

## 🚀 Usage

Clone this repository:

```bash
git clone https://github.com/LyricZL/Spintronics_Bristol.git
cd Spintronics_Bristol
```

Run simulations from MATLAB:

### 1. Simulate spin effects:
```matlab
cd simulation
simulate_spin_effects
```

### 2. Simulate spintronic devices:
```matlab
cd ../devices
simulate_devices
```

### 3. Compare integration methods:
```matlab
cd ../simulation
integration_comparison
```

## 📊 Sample Results

- Magnetisation vs. time in x/y/z under different materials and effects
- Performance comparison: RK4 > ode45 > AB2 (tradeoff between precision and speed)
- Device response analysis: STNO shows high-frequency oscillation, MTJ offers stability, Spin Valve shows isotropic behavior

## 📎 Notes

- All physical constants and device parameters are based on experimental literature and material databases.
- Figures are saved in `/plots` after execution.
- You can easily adjust simulation duration, material selection, and THz field parameters in the script headers.

## 📖 Background

This work investigates the application of spintronic materials and devices in terahertz communication systems, focusing on magnetisation dynamics, spin wave propagation, and response analysis under different effects and frequencies.

## 📝 License

This project is open-sourced under the MIT License. See the [LICENSE](LICENSE) file for details.

---
