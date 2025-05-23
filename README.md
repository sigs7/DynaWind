# DynaWind

**DynaWind** is a modular, high-fidelity Python simulation framework for co-simulation of wind turbine dynamics and power systems. It integrates aerodynamic and structural models from OpenFAST (FMU-based) with detailed electrical models of the turbine generator and converters.

## Features

- ⚙️ **Co-simulation with OpenFAST**: FMU integration for high-fidelity turbine dynamics.
- ⚡ **Permanent Magnet Synchronous Machine (PMSM)**: Detailed electric drive modeling including current control.
- 🔋 **DC-link System**: Voltage regulation with chopper logic and anti-windup control.
- 🔌 **Grid-Side Converter (PQ/PV control)**: Integrated with external power system solvers.
- 📈 **Logging and Visualization**: Automated results collection and flexible plotting with Matplotlib and Plotly.

## Directory Structure

```
cosim_models/
│
├── controller.py        # PI(D) controllers with anti-windup
├── dclink.py            # DC-link voltage dynamics and chopper logic
├── fast.py              # FMU handling for OpenFAST
├── ideal_generator.py   # Optional simplified generator model
├── pi_controller.py     # Alternative implementation of PI control
├── pmsm.py              # Permanent Magnet Synchronous Machine model
├── results.py           # Logging and plotting of simulation data
├── windturbine.py       # Top-level wind turbine model integrating all components
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
- `matplotlib`
- `plotly`
- `fmpy`
- `pandas`

