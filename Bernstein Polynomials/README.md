Bernstein Polynomials
=====================

I first saw the Bernstein polynomials used in a proof of Weierstrass' approximation theorem. Despite their relatively slow convergence rate for approximating functions, these polynomials have several interesting properties:
-The approximation for continuous functions is uniform, not pointwise.
-Derivatives of Bernstein polynomials converge to the approximated function's derivatives.
-The sign of the derivatives agrees with the approximated function's derivatives, so the Bernstein expansion is concave, convex, increasing, decreasing in the same regions as the original function.
-De Casteljau's algorithm allows for numerically stable evaluation of the Bernstein polynomial expansions.
