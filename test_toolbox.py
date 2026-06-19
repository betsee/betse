import numpy as np
import torch
from betse.science import sim_toolbox as stb
class MockParams:
    F = 96485.33  # Faraday constant
    R = 8.314     # Gas constant

print("🧪 KSL Sovereign Toolbox Verification")
print("-------------------------------------")

# Generate massive dummy data (1 Million points)
N = 1_000_000
print(f"   -> Generating {N} data points...")
cA = np.random.rand(N).astype(np.float32)
cB = np.random.rand(N).astype(np.float32)
vBA = np.random.rand(N).astype(np.float32) * 0.1 # Small voltages
Dc = np.ones(N).astype(np.float32) * 1e-9
d = np.ones(N).astype(np.float32) * 1e-6
zc = 1.0
T = 300.0
p = MockParams()

print("⚡ IGNITION: Calling electroflux (MPS)...")
try:
    flux = stb.electroflux(cA, cB, Dc, d, zc, vBA, T, p)
    print(f"✅ Success! Flux calculated.")
    print(f"   Shape: {flux.shape}")
    print(f"   Sample Value: {flux[0]}")
except Exception as e:
    print(f"❌ Failed: {e}")

