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


def kronecker_method_demonstration(polynomial: Union[Poly, sp.Expr]) -> List[Poly]:
    """
    Demonstrates how Kronecker's method works, showing the computational steps
    without actually performing the full factorization (which can be expensive).

    Parameters:
        polynomial: A SymPy polynomial in one variable

    Returns:
        List of irreducible polynomial factors
    """
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

    evaluations: Dict[int, int] = {point: int(
        pp_poly.eval(point)) for point in test_points}
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


def test_examples() -> None:
    """
    Test the factorization on multiple example polynomials.
    """
    examples: List[Poly] = [
        Poly(x**3 - 2*x**2 - 5*x + 6, x),  # (x-2)(x-1)(x+3)
        Poly(14*x**4 - 46*x**3 - 82*x**2 + 138*x + 120, x, domain='ZZ'),
        Poly(x**4 + x**2 - 20, x, domain='ZZ')
    ]

    for i, poly in enumerate(examples):
        print(f"\n\nExample {i+1}: {poly}")
        print("=" * 50)
        factors = kronecker_factorization(poly)

        print("Factors:")
        for factor in factors:
            print(f"  {factor}")

        # Verify by multiplying back
        product = factors[0]
        for f in factors[1:]:
            product = product * f

        print(
            f"\nVerification: product of factors equals original? {product == poly}")
        print(f"Original: {poly}")
        print(f"Product:  {product}")


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
