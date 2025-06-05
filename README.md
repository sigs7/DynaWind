# DynaWind

**DynaWind** is a modular, high-fidelity Python framework for simulating the dynamic interactions between modern wind turbines and electrical power systems. The framework combines detailed electrical modeling (generator, converter, DC-link) with high-fidelity structural and aerodynamic simulation through OpenFAST FMU co-simulation.

---

## 🔧 Features

- ⚙️ **Co-simulation with OpenFAST**  
  Integrates turbine aerodynamics and structural response via FMU (Functional Mock-up Unit).

- ⚡ **Permanent Magnet Synchronous Machine (PMSM)**  
  Includes dq0-based dynamic modeling and current control.

- 🔋 **DC-link System**  
  Supports voltage regulation, chopper logic, and anti-windup protection.

- 🔌 **Grid-Side Converter (GSC)**  
  Supports PQ and PV control modes; compatible with RMS-based power system models.

- 📊 **Logging and Visualization**  
  Simulation results are stored and visualized using Matplotlib and Plotly.

---

## 📁 Directory Structure

```
DynaWind/
└── dynawind/
    ├── dynawind_models/
    │   ├── __pycache__/
    │   ├── __init__.py
    │   ├── controller.py         # PI control with anti-windup
    │   ├── dclink.py             # DC-link voltage dynamics and chopper logic
    │   ├── fast.py               # FMU wrapper for OpenFAST
    │   ├── ideal_generator.py    # Optional simplified generator
    │   ├── pmsm.py               # PMSM modeling and control
    │   ├── results.py            # Logging and results export
    │   └── windturbine.py        # Top-level wind turbine integration
    │
    ├── figures/
    │   ├── Paper_results_60_SC/
    │   ├── Paper_results_120/
    │   └── Paper_results_360/
    │
    ├── plotting.py               # Custom plotting routines
    └── simulation.py             # Example simulation setup
```

---

## 🚀 Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/laherman1908/DynaWind.git
cd DynaWind
```

### 2. Install Dependencies

Create a virtual environment (recommended):

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

Install required packages:

```bash
pip install -r requirements.txt
```

### 3. Configure FMU

Download or generate an [OpenFAST](https://github.com/OpenFAST/openfast) FMU compatible with your wind turbine case.  
Update the FMU path in `fast.py` or through the configuration section in `windturbine.py`.

---

## 📦 Dependencies

- `numpy`
- `scipy`
- `fmpy`
- `pandas`
- `matplotlib`
- `plotly`
- `ipympl`


---

## 🧪 Run a Simulation

To run a predefined example simulation:

```bash
python simulation.py
```

Output data and figures will be saved to the `figures/` directory.

---

## 🤝 Acknowledgments

This work was developed as part of a master’s thesis at NTNU Trondheim, integrating OpenFAST and TOPS to study the dynamic behavior of wind turbines in future power systems.

---

## 🔗 Links

- [OpenFAST Repository](https://github.com/OpenFAST/openfast)  
- [FMU Standard - FMI](https://fmi-standard.org/)
- [TOPS Repository](https://github.com/hallvar-h/TOPS)
