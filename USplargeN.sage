# Some analytic expressions, valid at large N, for the number of
# singlets in the k-th tensor power of the n-th Adams operator applied
# to the fundamental of USp.

# Double factorial (of an odd number).
def df(p):
    assert (p%2 == 1)
    t = (p+1)/2
    return factorial(2*t)/(2^t*factorial(t))

def __I_n_even(n, k):
    result = 0
    for p in range(k/2+1):
        result += binomial(k, 2*p)* df(2*p-1) * n^p
    return (-1)^k*result

def __I_n_odd(n, k):
    if k%2 == 1:
        return 0
    return df(k-1)*n^(k/2)
    
def I(n, k):
    # Even and odd cases are different enough that it is better to
    # treat them separately.
    if n % 2 == 0:
        return __I_n_even(n, k)
    return __I_n_odd(n, k)

# The large N expression for the number of singlets in a k-th power of
# the n-th Adams operator applied to the adjoint.
def SAdj(n, k):
    result = 0
    for p in range(k+1):
        result += binomial(k, p) * I(n,2*p) * I(2*n,(k-p))
    return result/(2^k)
