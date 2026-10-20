#!/usr/bin/env python3
"""Golden vectors for fp8_mac: exhaustive aligned products and random dot products."""
import random
from fractions import Fraction

def e4m3(v):
    s, e, m = v >> 7, (v >> 3) & 15, v & 7
    val = Fraction(m, 8) * Fraction(2) ** (-6) if e == 0 else (1 + Fraction(m, 8)) * Fraction(2) ** (e - 7)
    return -val if s else val

UNIT = Fraction(1, 2 ** 18)
with open("build/fp8_products.hex", "w") as f:
    for k in range(65536):
        w, x = k & 255, k >> 8
        q = e4m3(w) * e4m3(x) / UNIT
        assert q.denominator == 1
        f.write(f"{int(q) & ((1 << 37) - 1):010x}\n")
random.seed(7)
with open("build/fp8_dots.hex", "w") as f:
    for t in range(200):
        L = random.choice([1, 16, 64, 256, 1024])
        ws = [random.randrange(256) for _ in range(L)]
        xs = [random.randrange(256) for _ in range(L)]
        tot = sum(e4m3(w) * e4m3(x) for w, x in zip(ws, xs)) / UNIT
        f.write(f"{L:04x}\n")
        for w, x in zip(ws, xs):
            f.write(f"{w:02x}{x:02x}\n")
        f.write(f"{int(tot) & ((1 << 48) - 1):012x}\n")
print("vectors written")
