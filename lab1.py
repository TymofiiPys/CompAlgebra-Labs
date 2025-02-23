import sympy as sp
from sympy.abc import x
from sympy.polys import Poly
from sympy import div, ZZ
import math
from itertools import product
from typing import Tuple, Union, List, Set


def divisors(n: Union[int, sp.Integer, float]) -> Set[int]:
    """
    Пошук всіх дільників числа n.
    """
    if n == 0:
        return [0]

    # Конвертація в ціле число
    if hasattr(n, 'evalf'):
        n = int(float(n))
    else:
        n = int(n)

    n = abs(n)
    divisors = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divisors.add(i)
            divisors.add(-i)
            if i != n // i:
                divisors.add(n // i)
                divisors.add(- (n // i))

    return divisors

def interpol()

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
        dy.append(divisors(y[i]))

    m = list(product(*dy))
    
    for mu in m:
        g = interpol()