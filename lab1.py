import sympy as sp
from sympy.abc import x
from sympy.polys import Poly
from sympy import div, ZZ
import math
from typing import Tuple

def divisors

def kronecker(f : Poly) -> Tuple[int, Poly, bool, Poly, Poly]:
    n = f.degree()
    if n == 0:
        return sp.Abs(f), sp.sign(f), False, None, None
    
    cf, prf = f.primitive()

    m = math.floor(n/2)

    x = range(m + 1)
    y = []
    dy = []

    for i in range(m + 1):
        y.append(f.subs(x, i))

        if y[i] == 0:
            return cf, prf, False, x - x[i], div(prf, x - x[i], domain=ZZ)
        dy.append()
        