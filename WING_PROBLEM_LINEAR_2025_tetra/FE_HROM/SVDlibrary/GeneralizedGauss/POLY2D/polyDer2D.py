import numpy as np

def polyDer2D(p, x, y, n, m):
    """
    polyDer2D(p, x, y, n, m)

    Evaluate derivatives of a 2-D polynomial using Horner's method.

    Evaluates the derivatives of 2-D polynomial `p` at the points
    specified by `x` and `y`, which must be the same dimensions. The
    outputs `(fx, fy)` will have the same dimensions as `x` and `y`.
    The order of `x` and `y` are specified by `n` and `m`, respectively.

    Parameters
    ----------
    p : array_like
        Polynomial coefficients in order specified by polyVal2D.html.
    x : array_like
        Values of 1st independent variable.
    y : array_like
        Values of 2nd independent variable.
    n : int
        Order of the 1st independent variable, `x`.
    m : int
        Order of the 2nd independent variable, `y`.

    Returns
    -------
    fx : ndarray
        Derivative with respect to x.
    fy : ndarray
        Derivative with respect to y.

    See Also
    --------
    numpy.polynomial.polynomial.polyval2d : Evaluate a 2-D polynomial at points (x, y).

    Example
    --------
    >>> print polyVal2D([1,2,3,4,5,6],2,3,2,1)
    >>> (39, 11)
    >>> 1*2*2*3 + 2*3 + 4*2*2 + 5
    39
    >>> 1*(2**2) + 2*2 + 3
    11
    >>> print polyVal2D([1,2,3,4,5,6,7,8,9],2,3,2,2)
    >>> (153, 98)
    >>> 1*2*2*(3**2) + 2*(3**2) + 4*2*2*3 + 5*3 + 7*2*2 + 8
    153
    >>> 1*2*(2**2)*3 + 2*2*2*3 + 3*2*3 + 4*(2**2) + 5*2 + 6
    98
    """
    # TODO: check input args
    p = np.array(p)
    x = np.array(x)
    y = np.array(y)
    n = np.array(n)
    m = np.array(m)

    # fx = df/dx
    fx = n * p[0]
    for ni in np.arange(n - 1):
        fx = fx * x + (n - ni - 1) * p[1 + ni]
    for mi in np.arange(m):
        mj = (n + 1) * (mi + 1)
        gx = n * p[mj]
        for ni in np.arange(n - 1):
            gx = gx * x + (n - ni - 1) * p[mj + 1 + ni]
        fx = fx * y + gx

    # fy = df/dy
    fy = p[0]
    for ni in np.arange(n):
        fy = fy * x + p[1 + ni]
    fy = m * fy
    for mi in np.arange(m - 1):
        mj = (n + 1) * (mi + 1)
        gy = p[mj]
        for ni in np.arange(n):
            gy = gy * x + p[mj + 1 + ni]
        fy = fy * y + (m - mi - 1) * gy

    return fx, fy
