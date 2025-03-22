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


def interpol(xl: List[int], yl: List[int]):
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


def kronecker(f: Poly) -> Tuple[int, Poly, bool, Poly, Poly]:
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
            return cf, prf, False, Poly(x - xl[i]), div(prf, x - xl[i])[0]
        dy.append(divisors(yl[i]))

    m = list(product(*dy))

    for mu in m:
        g = interpol(list(xl), mu)

        if f.domain == ZZ and g.degree() >= 1 and rem(prf, g) == 0 and div(prf, g)[0].domain == ZZ:
            if LC(g) < 0:
                g = -g
            return cf, prf, False, g, div(prf, g)[0]

    return cf, prf, True, prf, 1


def kroneckerIRR(f: Poly):
    kr = kronecker(f)
    if kr[2]:
        return kr
    else:
        cf = kr[0]
        prf = kr[1]
        h = kr[3]
        while not kronecker(h)[2]:
            h = kronecker(h, x)[3]
        return cf, prf, False, h, div(prf, h)[0]


def kroneckerFACTS(f: Poly):
    n = f.degree()
    if n <= 1:
        return f
    
    scal, prf = f.primitive()
    kr = kroneckerIRR(prf)
    g, h = kr[3], kr[4]
    irr, mlt = [g], [1]
    
    while not (h == 1):
        kr = kroneckerIRR(h)
        g, h = kr[3], kr[4]
        try:
            k = irr.index(g)
            mlt[k] = mlt[k] + 1
        except ValueError:
            irr.append(g)
            mlt.append(1)
    
    lng = len(irr)
    ans = []

    for i in range(lng):
        ans.append((irr[i], mlt[i]))

    return scal, ans



if __name__ == "__main__":
    polynomials = [Poly(14*x**4 - 46*x**3 - 82*x**2 + 138*x + 120),
                   Poly(x**4 + x**2 - 20)]
    for polynomial in polynomials:
        print()
        print("Поліном:", polynomial)
        factors = kroneckerFACTS(polynomial)
        print("Скалярний множник:", factors[0])
        print("Незвідні множники:")
        for factor in factors[1]:
            print(factor[0], ", кратність:", factor[1])
        product_check = factors[0]
        for f in factors[1]:
            product_check = product_check * f[0] ** f[1]
        print("Перевірка множників: Добуток отриманих множників:", product_check)
        print("\n===================")
