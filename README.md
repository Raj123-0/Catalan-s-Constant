[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

===============================================================================
PROJECT: Catalan's Constant Computation Engine
===============================================================================

OVERVIEW:
Calculates Catalan's constant (G ≈ 0.91596559417721901505...) to arbitrary 
precision (N digits) using a fast-converging Ramanujan-Guillera hypergeometric 
series, parallel binary splitting, multi-core multiprocessing, and C-accelerated 
gmpy2 arithmetic.

ALGORITHM & MATHEMATICS:
- Ramanujan-Guillera Series:
    G = (1/8) * sum_{k=0}^{infinity} (-1)^k * (2k)!^3 / ((4k)! * 2^(4k)) * (18k^2 + 8k + 1) / (2k + 1)^2
- Fast Hypergeometric Convergence: ~2.408 decimal digits per evaluated series term.
- Parallel Binary-Splitting: Evaluates integer products across CPU worker pools.

## Usage

```bash
python "Catalan's Constant.py" --help
```
