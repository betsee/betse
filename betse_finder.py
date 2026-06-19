import os
import inspect
from betse.science.sim.simrunner import SimRunner
# Try to import the likely math/solver locations
try:
    import betse.science.math as math_pkg
    print(f"📍 Math Package Location: {os.path.dirname(math_pkg.__file__)}")
except ImportError:
    print("⚠️ Could not import betse.science.math directly")

print(f"📍 SimRunner Location: {inspect.getfile(SimRunner)}")

# Let's list the files in the science/math directory if we can find it
try:
    target_dir = os.path.join(os.path.dirname(inspect.getfile(SimRunner)), '../math')
    target_dir = os.path.normpath(target_dir)
    print(f"\n📂 Scanning: {target_dir}")
    files = [f for f in os.listdir(target_dir) if f.endswith('.py')]
    print("   Found:", files)
except Exception as e:
    print(f"   Scan failed: {e}")

