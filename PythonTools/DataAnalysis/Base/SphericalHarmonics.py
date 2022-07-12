
import math

# Calculate combinatorics
def n_choose_r(n : int, r : int):
    assert (r >= 0) and (n >= r)
    return math.factorial(n) / (math.factorial(r) * math.factorial(n - r))


def spin_Y_lm(es : int, el : int, em : int, x : float = None, y : float = None, z : float = None,
              theta : float = None, phi : float = None):
    assert (el >= 0) and (el >= abs(em))

    # calculate useful position quantities
    if x is None:
        assert y is None and z is None
        assert theta is not None and phi is not None
    if theta is None:
        assert phi is None
        assert x is not None and y is not None and z is not None

        r = max(math.sqrt(x * x + y * y + z * z), 1e-6)
        theta = math.acos(z / r)
        phi = math.atan2(y, x)

    coefficient = math.pow(-1.0, es) * math.sqrt((2.0 * el + 1.0) / (4.0 * math.pi));
    coefficient *= math.sqrt(math.factorial(el + em) * math.factorial(el - em) /
                        math.factorial(el + es) / math.factorial(el - es))

    sum_terms = 0.0
    lower_limit = em + es if (em + es > 0) else 0
    upper_limit = el + em if (el + em < el + es) else el + es

    for i in range(lower_limit, upper_limit+1):
        temp = n_choose_r(el + es, i) * n_choose_r(el - es, i - es - em)
        sum_terms += temp * math.pow(-1.0, i) *\
               math.pow(math.cos(theta / 2.0), 2 * (el - i) + es + em) *\
               math.pow(math.sin(theta / 2.0), 2 * i - em - es)

    y_lm_Real = coefficient * sum_terms * math.cos(em * phi)
    y_lm_Im = coefficient * sum_terms * math.sin(em * phi)

    return y_lm_Real + 1j * y_lm_Im

def check_integral(es1, el1, em1, es2, el2, em2):
    ntheta = 50
    nphi = 50
    total = 0
    dtheta = math.pi / (ntheta-1)
    dphi = 2 * math.pi / nphi
    for t in range(ntheta):
        theta = t * dtheta
        for p in range(nphi):
            phi = p * dphi
            area_alement = math.sin(theta) * dtheta * dphi
            total += area_alement * spin_Y_lm(theta=theta, phi=phi, es=es1, el=el1, em=em1) * spin_Y_lm(theta=theta, phi=phi, es=es2, el=el2, em=em2).conjugate()
    return total

# run this file to just run some checks
if __name__ == "__main__":
    print(spin_Y_lm(theta=1, phi=1, es=-2, em=2, el=2))
    print(spin_Y_lm(theta=1, phi=1, es=0, em=2, el=2))

    integral = check_integral(1,1,1, 1,1,1)
    print(integral.real)
    print(integral.imag)
