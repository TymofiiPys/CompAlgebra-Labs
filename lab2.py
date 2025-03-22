import sympy as sp
from sympy.abc import x, y, z
from sympy.polys import Poly
from typing import Set

# Варіант 3

def reduced_groebner_basis(F : Set[Poly]):
    G = list(groebner_basis(F))
    G_monic = [sp.monic(gi) for gi in G]
    G_minimal = set()
    removed_ind = []
    for i in range(len(G_monic)):
        gi = G_monic[i]
        print("gi: ", gi)
        lmgi = sp.LM(gi)
        divisible_by_any = False
        for j in range(len(G_monic)):
            if j != i and j not in removed_ind:
                lmgj = sp.LM (G_monic[j])
                if lmgi % lmgj == 0:
                    removed_ind.append(i)
                    divisible_by_any = True
                    break
        if not divisible_by_any:
            G_minimal.add(gi)
    
    print(G_minimal)
    print()
    G_reduced = []

    for gi in G_minimal:
        _, gri = sp.reduced(gi, G_minimal.difference({gi}))
        if gri != 0:
            G_reduced.append(gri)
    
    return G_reduced
    

def spoly(f : Poly, g : Poly) -> Poly:
    ltf = sp.LT(f)
    ltg = sp.LT(g)
    return sp.lcm(ltf, ltg) * (f/ltf - g/ltg)

def groebner_basis(F : Set[Poly]) -> Set[Poly]:
    G = F
    C = {frozenset([g1, g2]) for g1 in G for g2 in G if not g1 == g2} 
    while len(C) > 0:
        f, g = C.pop()
        print("f and g:")
        print(f, g)
        _, h = sp.reduced(spoly(f, g), G)
        print("h:")
        print(h)
        print()
        if not h == 0:
            C.update({frozenset([g1, h]) for g1 in G if not g1 == h})
            G.add(h)
    return G

if __name__ == "__main__":
    system = {Poly(x*z-2*y+1),
              Poly(y*z - 1 + z),
              Poly(y*z+x*y*z+z)}
    G = reduced_groebner_basis(system)
    print(G)
    print(sp.groebner(system, x, y, z, order='grlex'))
    
