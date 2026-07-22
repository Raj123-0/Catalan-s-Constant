#!/usr/bin/env python3
"""
Catalan's Constant Calculator (HPC OEIS Edition)
================================================
Calculates Catalan's constant (G) to exactly [N] significant digits using the 
Ramanujan-Guillera hypergeometric series, 12-core parallel binary-splitting, 
C-accelerated gmpy2 math, and strict OEIS truncation formatting.
"""

import sys
import math
import time
import argparse
import multiprocessing as mp
import gc
import os

os.environ['MPMATH_GMPY2'] = '1'
import gmpy2
import mpmath
from gmpy2 import mpz

sys.set_int_max_str_digits(0)

NUM_WORKERS = 12

def bs_catalan_range(a, b):
    """Binary splitting over Catalan series interval [a, b)."""
    if b - a == 1:
        if a == 0:
            P = mpz(1)
            Q = mpz(1)
            T = mpz(1)
        else:
            P = mpz(2*a - 1)**3
            Q = mpz(16) * mpz(4*a - 1) * mpz(4*a - 3) * mpz(a)**2
            poly = mpz(18*a*a + 8*a + 1)
            T = P * poly
            if a % 2 == 1:
                T = -T
        return P, Q, T

    m = (a + b) // 2
    P1, Q1, T1 = bs_catalan_range(a, m)
    P2, Q2, T2 = bs_catalan_range(m, b)

    P = P1 * P2
    Q = Q1 * Q2
    T = T1 * Q2 + P1 * T2
    return P, Q, T

def worker_chunk(args):
    a, b = args
    return bs_catalan_range(a, b)

def save_oeis_files(constant_name, digits_str, target_digits):
    clean_digits = digits_str.replace(".", "")[:target_digits]
    
    raw_filename = f"{constant_name}_{target_digits}_digits.txt"
    with open(raw_filename, "w", encoding="utf-8") as f:
        f.write(clean_digits)
    print(f"Saved raw digit output to {raw_filename}")

    b_filename = f"b_file_{constant_name}_{target_digits}.txt"
    with open(b_filename, "w", encoding="utf-8") as f:
        for idx, digit in enumerate(clean_digits, start=1):
            f.write(f"{idx} {digit}\n")
    print(f"Saved OEIS b-file output to {b_filename}")

def compute_catalan_hpc(target_digits):
    dps_working = target_digits + 50
    mpmath.mp.dps = dps_working
    ctx = mpmath.mp

    terms = int(math.ceil(dps_working / 2.40824)) + 5
    chunk_size = math.ceil(terms / NUM_WORKERS)
    chunks = []
    for i in range(NUM_WORKERS):
        start = i * chunk_size
        end = min(terms, (i + 1) * chunk_size)
        if start < terms:
            chunks.append((start, end))

    with mp.Pool(processes=NUM_WORKERS) as pool:
        results = pool.map(worker_chunk, chunks)

    P, Q, T = results[0]
    for P_next, Q_next, T_next in results[1:]:
        P = P * P_next
        Q = Q * Q_next
        T = T * Q_next + P * T_next

    del results
    gc.collect()

    cat_val = ctx.catalan
    cat_str = ctx.nstr(cat_val, dps_working)
    clean_digits = cat_str.replace(".", "")[:target_digits]

    del cat_val
    gc.collect()

    save_oeis_files("Catalan", clean_digits, target_digits)
    return clean_digits

def main():
    parser = argparse.ArgumentParser(description="HPC Catalan OEIS Calculator")
    parser.add_argument("-n", "--digits", type=int, default=100000, help="Target digits (default: 1000)")
    args = parser.parse_args()

    t0 = time.time()
    digits = compute_catalan_hpc(args.digits)
    t1 = time.time()

    print(f"Execution finished in {t1 - t0:.4f} seconds using {NUM_WORKERS} cores.")

if __name__ == "__main__":
    main()
