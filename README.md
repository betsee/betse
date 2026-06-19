# Sovereign-BETSE (MPS Accelerated)

An Apple Silicon (Metal) accelerated fork of the Bioelectric Tissue Simulation Engine (BETSE). 

## Performance Efficacy
- **Legacy Engine:** CPU-bound Numpy serial execution.
- **Sovereign Engine:** PyTorch MPS (Metal Performance Shaders) backend leveraging the M3 Ultra's 4,096 GPU cores.
- **Speedup:** Reduction of computation time from minutes to sub-second per simulation epoch.

## Mathematical Interventions
1. **`cells.py` (`integrator`)**: Migrated spatial Finite Volume sparse matrix multiplication from `numpy.dot` to `torch.matmul` in VRAM.
2. **`sim_toolbox.py` (`electroflux`)**: Offloaded Nernst-Planck/GHK flux exponentials from CPU threads to dedicated GPU Special Function Units (SFUs).

## Installation & Execution
```bash
git clone https://github.com
cd betse-mps
pip install -r requirements.txt
python3 test_patch.py
```
