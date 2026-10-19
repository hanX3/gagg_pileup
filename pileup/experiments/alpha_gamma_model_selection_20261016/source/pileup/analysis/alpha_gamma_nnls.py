"""Equivalent tall-matrix NNLS reduction, with numerical safety checks.

R.T R = A.T A and R.T z = A.T b imply that ||R x-z||² and
||A x-b||² differ by a constant. No penalty, prior or threshold is added.
Small or numerically unsuitable systems use SciPy's original solver.
"""
import numpy as np
from scipy.linalg import cholesky, solve_triangular, qr, LinAlgError
from scipy.optimize import nnls as direct_nnls


def _qr_fallback(a, b, maxiter):
    """Avoid squaring the condition number when a tall Gram factor is unsafe."""
    try:
        orthogonal, upper = qr(a, mode='economic', check_finite=False)
        x, _ = direct_nnls(upper, orthogonal.T @ b, maxiter=maxiter)
        residual = a @ x-b
        gradient = a.T @ residual
        scale = max(1., float(np.max(np.abs(a.T @ b))))
        active = x > 0
        if (np.min(gradient) >= -1e-8*scale and
                (not active.any() or np.max(np.abs(gradient[active])) <= 1e-8*scale)):
            return x, float(np.linalg.norm(residual))
    except (LinAlgError, RuntimeError):
        pass
    return direct_nnls(a, b, maxiter=maxiter)


def reduced_nnls(a, b, *, maxiter=None):
    if a.shape[1] < 16 or a.shape[0] < 4*a.shape[1]:
        return direct_nnls(a, b, maxiter=maxiter)
    gram = a.T @ a
    rhs = a.T @ b
    try:
        upper = cholesky(gram, lower=False, check_finite=False)
        # A nearly rank-deficient reduction can amplify rounding error.
        if np.min(np.diag(upper)) < 1e-6*np.max(np.diag(upper)):
            return _qr_fallback(a, b, maxiter)
        z = solve_triangular(upper.T, rhs, lower=True, check_finite=False)
        x, _ = direct_nnls(upper, z, maxiter=maxiter)
        residual = a @ x-b
        gradient = a.T @ residual
        scale = max(1., float(np.max(np.abs(rhs))))
        active = x > 0
        if (np.min(gradient) < -1e-8*scale or
                (active.any() and np.max(np.abs(gradient[active])) > 1e-8*scale)):
            return _qr_fallback(a, b, maxiter)
        return x, float(np.linalg.norm(residual))
    except (LinAlgError, RuntimeError):
        return _qr_fallback(a, b, maxiter)
