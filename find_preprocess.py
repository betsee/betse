import os

print("🔍 Scanning local BETSE codebase for 'class Preprocess'...")
# Walk through the directory to find the file defining the class
for root, dirs, files in os.walk("."):
    for file in files:
        if file.endswith(".py"):
            path = os.path.join(root, file)
            try:
                with open(path, "r", encoding="utf-8") as f:
                    if "class Preprocess" in f.read():
                        # Convert file path to Python import path
                        clean_path = path.replace("./", "").replace("/", ".").replace(".py", "")
                        print(f"\n✅ FOUND IT! File: {path}")
                        print(f"   Likely Import Path: {clean_path}")
            except Exception:
                pass
print("\nScan complete.")
