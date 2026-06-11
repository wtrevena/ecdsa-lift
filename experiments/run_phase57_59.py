import runpy, sys, pathlib, traceback
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
for mod in ["experiments.phase57_pseudovariance",
            "experiments.phase58_carry_lemma",
            "experiments.phase59_u3"]:
    print("\n========== RUN", mod, "==========")
    try:
        runpy.run_module(mod, run_name="__main__")
    except Exception:
        traceback.print_exc()
print("\nALL DONE")
