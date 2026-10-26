import pickle, sys, time
from .data import generate
for name in sys.argv[1:]:
    t = time.time()
    g, d = generate(name, T=6000, seed={"ieee14": 1, "ieee118": 2}[name])
    pickle.dump(d, open(f"data/{name}.pkl", "wb"))
    print(name, "done", time.time() - t, flush=True)
