import sympy as sp
from sympy.abc import x
from sympy.polys import Poly
from sympy import div, rem, ZZ, simplify, LC
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

def interpol(xl : List[int], yl : List[int]):
    n = len(xl)
    if n == 0 or n != len(yl) or not len(xl) == len(set(xl)):
        print("Error")
        exit(1)

    f = Poly(0, x)
    for i in range(n):
        p = yl[i]
        for j in range(n):
            if j != i:
                p = simplify(p * (x - xl[j]) / (xl[i] - xl[j]))
        
        f = simplify(f + p)

    return f


def kronecker(f : Poly) -> Tuple[int, Poly, bool, Poly, Poly]:
    n = f.degree()
    if n == 0:
        return sp.Abs(f), sp.sign(f), False, None, None
    
    cf, prf = f.primitive()

    m = math.floor(n/2)

    xl = range(m + 1)
    yl = []
    dy = []

    for i in range(m + 1):
        yl.append(f.subs(x, xl[i]).evalf())

        if yl[i] == 0:
            return cf, prf, False, x - xl[i], div(prf, x - xl[i])[0]
        dy.append(divisors(yl[i]))

    m = list(product(*dy))
    
    for mu in m:
        g = interpol(list(xl), mu)

        if f.domain == ZZ and g.degree() >= 1 and rem(prf, g) == 0 and div(prf, g)[0].domain == ZZ:
            if LC(g) < 0:
                g = -g
            return cf, prf, False, g, div(prf, g)[0]
        
    return cf, prf, True, None, None

if __name__ == "__main__":
    #print(interpol([0,1,2,3], [-2,4,5,3]))        
    #print(kronecker(Poly(2*x**5 + x**4 + x**2 +x +2)))
    print(kronecker(Poly(18*x**5 + 6*x**4 + 12*x**3 + 27*x**2 + 9*x + 18)))
    #print(kronecker(Poly(2*x**4 + 8)))