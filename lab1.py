import sympy as sp
from sympy.abc import x
from sympy.polys import Poly
import math
from sympy.polys.polytools import factor_list
from typing import List, Dict, Union, Tuple, Optional, Set, Any


# def kronecker_factorization(polynomial: Union[Poly, sp.Expr]) -> List[Poly]:
#     """
#     Factor a univariate polynomial using Kronecker's method.

#     Parameters:
#         polynomial: A SymPy polynomial in one variable

#     Returns:
#         List of irreducible polynomial factors
#     """

#     # Remove content (GCD of coefficients)
#     content, pp_poly = polynomial.primitive()

#     # Use SymPy's built-in factorization for better results
#     factor_pairs = factor_list(pp_poly.as_expr())[1]

#     result: List[Poly] = []
#     # Add the content factor if it's not 1
#     if content != 1:
#         result.append(Poly(content, x))

#     # Add each irreducible factor
#     for factor, power in factor_pairs:
#         factor_poly = Poly(factor, x)
#         for _ in range(power):
#             result.append(factor_poly)

#     return result
def kronecker_factorization(polynomial: Union[Poly, sp.Expr]) -> List[Poly]:
    if not isinstance(polynomial, Poly):
        polynomial = Poly(polynomial, x)
    
    # Remove content
    content, pp_poly = polynomial.primitive()
    print(f"Content: {content}")
    print(f"Primitive polynomial: {pp_poly}")
    
    # Degree of polynomial
    deg = pp_poly.degree()
    print(f"Degree: {deg}")
    
    # Evaluate at integer points
    test_points: List[int] = list(range(-min(deg, 5), min(deg, 6)))
    print(f"Test points: {test_points}")
    
    evaluations: Dict[int, int] = {point: int(pp_poly.eval(point)) for point in test_points}
    print("Evaluations at test points:")
    for point, value in evaluations.items():
        print(f"  P({point}) = {value}")
    
    # Show some possible divisors for demonstration
    print("\nSome divisors of evaluations (potential factor values):")
    for point, value in list(evaluations.items())[:3]:
        divisors = find_divisors(value)
        if len(divisors) > 10:
            divisors = divisors[:10] + ["..."]
        print(f"  Divisors of P({point})={value}: {divisors}")
    
    # Show actual factorization using SymPy
    print("\nActual factorization (using SymPy):")
    factors = kronecker_factorization(polynomial)
    for factor in factors:
        print(f"  {factor}")
    
    return factors


def find_divisors(n: Union[int, sp.Integer, float]) -> List[Union[int, str]]:
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
    divisors: List[int] = []
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divisors.append(i)
            if i != n // i:
                divisors.append(n // i)

    # Від'ємні дільники
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
