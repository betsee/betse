import sys
import os
import numpy as np
import torch

print("🧪 KSL Sovereign Patch Verification Protocol")
print("------------------------------------------")

try:
    # 1. Import the class we patched
    from betse.science.cells import Cells
    print("✅ Imported Cells class successfully.")
except ImportError as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

# 2. Monkey-Patch the constructor (init) to skip the heavy setup
# We don't want to load config files, we just want to test the math engine.
def mock_init(self):
    print("   -> 🔧 Mocking Cells environment...")
    pass

Cells.__init__ = mock_init

# 3. Instantiate the subject
sim_cells = Cells()

# 4. Inject Dummy Data (Simulating a Grid)
# We need to mimic the data structures the integrator expects
num_cells = 100
num_mems = 400

print("   -> 🎲 Generating dummy physics data...")
# M_sum_mems: Connectivity matrix (Cells x Membranes)
sim_cells.M_sum_mems = np.random.rand(num_cells, num_mems).astype(np.float32)
# num_mems: Divisor
sim_cells.num_mems = np.ones(num_cells, dtype=np.float32) * 4.0
# Indices for neighbor lookup
sim_cells.nn_i = np.random.randint(0, num_mems, num_mems)
sim_cells.mem_i = np.arange(num_mems)

# 5. Prepare Inputs
f_center = np.random.rand(num_cells).astype(np.float32)
f_mem = np.random.rand(num_mems).astype(np.float32)

# 6. TRIGGER THE SOVEREIGN ENGINE
print("\n⚡ ATTEMPTING IGNITION (Calling integrator)...")
try:
    # This call should trigger your "KSL: Initializing..." print statement
    result_center, result_mem = sim_cells.integrator(f_center, f_mem)
    
    print("\n✅ Execution Successful!")
    print(f"   Output Shape: {result_center.shape}")
    
    # Check if a Torch tensor was secretly created in self
    if hasattr(sim_cells, '_M_sum_mems_mps'):
        print("   🏆 SOVEREIGN STATUS: CONFIRMED.")
        print(f"   Cache detected on device: {sim_cells._M_sum_mems_mps.device}")
    else:
        print("   ⚠️ WARNING: Code ran, but Sovereign Cache was not created.")

except Exception as e:
    print(f"\n❌ CRITICAL FAILURE in Integrator: {e}")
    import traceback
    traceback.print_exc()

print("------------------------------------------")

