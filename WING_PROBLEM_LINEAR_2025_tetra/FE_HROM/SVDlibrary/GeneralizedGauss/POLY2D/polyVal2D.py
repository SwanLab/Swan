import numpy as np

def polyVal2D(p, x, y, n, m):
    """
    polyVal2D(p, x, y, n, m)

    Evaluate a 2-D polynomial using Horner's method.

    Evaluates the 2-D polynomial `p` at the points specified by `x`
    and `y`, which must be the same dimensions. The output `f` will
    have the same dimensions as `x` and `y`. The order of `x` and `y`
    are specified by `n` and `m`, respectively.

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
    f : ndarray
        Values of the evaluated 2-D polynomial.

    See Also
    --------
    numpy.polynomial.polynomial.polyval2d : Evaluate a 2-D polynomial at points (x, y).

    Example
    --------
    >>> (5**2 * 6**2 + 2 * 5 * 6**2 + 3 * 6**2 + 4 * 5**2 * 6 +
    ...  5 * 5 * 6 + 6 * 6 + 7 * 5**2 + 8 * 5 + 9)
    2378
    >>> polyVal2D([1,2,3,4,5,6,7,8,9],5,6,2,2)
    2378
    >>> 5 * 6 + 2 * 6 + 3 * 5 + 4
    61
    >>> polyVal2D([1,2,3,4],5,6,1,1)
    61
    """
    # TODO: check input args
    p = np.array(p)
    x = np.array(x)
    y = np.array(y)
    n = np.array(n)
    m = np.array(m)
    
    # evaluate f(x,y)
    f = p[0]
    for ni in np.arange(n):
        f = f * x + p[1 + ni]
    for mi in np.arange(m):
        mj = (n + 1) * (mi + 1)
        g = p[mj]
        for ni in np.arange(n):
            g = g * x + p[mj + 1 + ni]
        f = f * y + g
    
    return f
