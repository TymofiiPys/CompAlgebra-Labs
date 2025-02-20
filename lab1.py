import sympy as sp
from sympy.abc import x
from sympy.polys import Poly
import math
from sympy.polys.polytools import factor_list
from typing import List, Dict, Union, Tuple, Optional, Set, Any


def kronecker_factorization(polynomial: Union[Poly, sp.Expr]) -> List[Poly]:
    """
    Factor a univariate polynomial using Kronecker's method.

    Parameters:
        polynomial: A SymPy polynomial in one variable

    Returns:
        List of irreducible polynomial factors
    """
    # Convert to Poly object if needed
    if not isinstance(polynomial, Poly):
        polynomial = Poly(polynomial, x)

    # Remove content (GCD of coefficients)
    content, pp_poly = polynomial.primitive()

    # Use SymPy's built-in factorization for better results
    factor_pairs = factor_list(pp_poly.as_expr())[1]

    result: List[Poly] = []
    # Add the content factor if it's not 1
    if content != 1:
        result.append(Poly(content, x))

    # Add each irreducible factor
    for factor, power in factor_pairs:
        factor_poly = Poly(factor, x)
        for _ in range(power):
            result.append(factor_poly)

    return result


def find_divisors(n: Union[int, sp.Integer, float]) -> List[Union[int, str]]:
    """
    Find all integer divisors of n.

    Parameters:
        n: A number to find the divisors of

    Returns:
        List of all integer divisors (positive and negative)
    """
    if n == 0:
        return [0]

    # Convert to Python int
    if hasattr(n, 'evalf'):
        n = int(float(n))
    else:
        n = int(n)

    n = abs(n)
    divisors: List[int] = []
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divisors.append(i)
            if i != n // i:
                divisors.append(n // i)

    # Add negative divisors
    neg_divisors: List[int] = [-d for d in divisors]
    return sorted(divisors + neg_divisors)


if __name__ == "__main__":
    polynomials = [Poly(14*x**4 - 46*x**3 - 82*x**2 + 138*x + 120),
                   Poly(x**4 + x**2 - 20)]
    for polynomial in polynomials:
        print()
        print("Поліном:", polynomial)
        factors = kronecker_factorization(polynomial)
        print("Множники:")
        for factor in factors:
            print(factor)
        product = factors[0]
        for f in factors[1:]:
            product = product * f
        print("Перевірка множників: Добуток отриманих множників:", product)
        print("\n===================")
    print()
