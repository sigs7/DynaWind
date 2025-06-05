# DynaWind

**DynaWind** is a modular, high-fidelity Python simulation framework for co-simulation of wind turbine dynamics and power systems. It integrates aerodynamic and structural models from OpenFAST (FMU-based) with detailed electrical models of the turbine generator and converters.

## Features

- ⚙️ **Co-simulation with OpenFAST**: FMU integration for high-fidelity turbine dynamics.
- ⚡ **Permanent Magnet Synchronous Machine (PMSM)**: Detailed electric drive modeling including current control.
- 🔋 **DC-link System**: Voltage regulation with chopper logic and anti-windup control.
- 🔌 **Grid-Side Converter (PQ/PV control)**: Integrated with external power system solvers.
- 📈 **Logging and Visualization**: Automated results collection and flexible plotting with Matplotlib and Plotly.


## 📁 Directory Structure
dynawind/
│
├── dynawind_models/ # Modular subcomponents of the wind turbine
│ ├── controller.py # Generic PI control schemes
│ ├── dclink.py # DC-link model and control
│ ├── fast.py # FMU wrapper for OpenFAST
│ ├── ideal_generator.py # Optional simplified generator
│ ├── pmsm.py # PMSM modeling and control
│ ├── results.py # Logging and results export
│ └── windturbine.py # System integrator model
│
├── figures/ # Output figures from simulations
│ ├── Paper_results_60_SC/
│ ├── Paper_results_120/
│ └── Paper_results_360/
│
├── plotting.py # Custom plotting routines
├── simulation.py # Example simulation setup
└── examples/ # Example scripts (optional)
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/laherman1908/DynaWind.git
  
   cd DynaWind
   ```

2. Install the dependencies (recommend using a virtual environment):
   ```bash
   pip install -r requirements.txt
   ```

3. Make sure you have a compatible OpenFAST FMU and update the FMU path in `windturbine.py`.

## Dependencies

- `numpy`
- `scipy`
- `fmpy`
- `matplotlib`
- `pandas`
- `ipympl`
- `plotly`

