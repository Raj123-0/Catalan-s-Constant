[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

===============================================================================
PROJECT: Catalan's Constant Computation Engine
===============================================================================

OVERVIEW:
Calculates Catalan's constant (G ≈ 0.91596559417721901505...) to arbitrary 
precision (N digits) using the Zucker-Bradley series, exact-integer binary 
splitting, and C-accelerated gmpy2 arithmetic.

ALGORITHM & MATHEMATICS:
- Zucker-Bradley Series:
    G = (pi/8) * ln(2 + sqrt(3)) + (3/8) * sum_{n>=0} 1 / ((2n+1)^2 * C(2n,n))
  where C(2n, n) is the central binomial coefficient.
- Geometric Convergence: terms decay like ~4^-n (about 1.7 terms per digit).
- Binary Splitting: the series is summed in exact integer arithmetic and
  collapses to a single fraction T/Q, requiring one final high-precision
  division.

NOTE: the previously documented "Ramanujan-Guillera" formula was
mathematically invalid (it diverges), and the previous code computed it
but discarded the result — digits came from mpmath's built-in catalan.
Both are now correct and actually used.

## Usage

```bash
python "Catalan's Constant.py" --help
```

TESTS:
    pytest tests/
