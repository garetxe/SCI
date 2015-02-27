# Some analytic expressions, valid at large N, for the number of
# singlets in the k-th tensor power of the n-th Adams operator applied
# to the fundamental of SO(2n+1). This works very similarly to the USp
# case, there is just some small sign differences.

# Double factorial (of an odd number).
def df(p):
    assert (p%2 == 1)
    t = (p+1)/2
    return factorial(2*t)/(2^t*factorial(t))

def __I_n_even(n, k):
    result = 0
    for p in range(k/2+1):
        result += binomial(k, 2*p)* df(2*p-1) * n^p
    return result

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
        result += (-1)^(k-p) * binomial(k, p) * I(n,2*p) * I(2*n,(k-p))
    return result/(2^k)

# The letter for SO(\infty) of degree n (not including the 1/n factor).
def i_n(n, t, x, xi, gs, rs, ps, pis, Fs, Fcs):
    # Gamma_n and gamma_{2n}
    gn = gs[n-1]
    g2n = gs[2*n-1]
    # Modulus and phase of rho_n and rho_{2n}
    rn = rs[n-1]
    pn = ps[n-1]
    pin = pis[n-1]
    r2n = rs[2*n-1]
    p2n = ps[2*n-1]
    pi2n = pis[2*n-1]
    # The flavor character
    Xn = Fs[n-1]
    Xcn = Fcs[n-1]

    # Vector multiplet
    result = (2*t^2 - t*(x + xi))*(1/2*(gn^2 - g2n) + rn^2 - 1)
    # Bifundamental
    result += t^(2/3)*gn*rn*pin*Xn - t^(4/3)*gn*rn*pn*Xcn
    # Antisymmetric
    result += t^(2/3)*(1/2)*((rn*pn)^2 - (r2n*p2n)) * Xn - \
        t^(4/3)*(1/2)*((rn*pin)^2 - (r2n*pi2n))*Xcn

    return result/((1-t*x)*(1-t*xi))

# Given a symbolic expression in g_n, r_n, do the integration.
def largeN_integration(m, maxn):
    # Treat the constant monomial separately, since sage does not
    # treat it as an element of a polynomial ring.
    if m in QQ:
        return (1, m)

    R = m.parent()
    names = R.variable_names()

    # The remaining monomial, after integrating over the gauge
    # group. It will be a funcion of the flavor characters only.
    mon = 1

    # The result of integrating over the gauge group
    c = 1 

    # Keep track of the degrees in the phases
    p_degs = [0]*int(2*maxn)

    for i, degree in enumerate(m.degrees()):
        if degree != 0:
            name = names[i]
            if name.startswith("g"):
                n = int(name[1:])
                c *= I(n, degree)
            elif name.startswith("r"):
                n = int(name[1:])
                k = degree/2
                if k in ZZ:
                    c *= factorial(k) * (n^k)
                else:
                    c = 0
            elif name.startswith("pi"):
                n = int(name[2:])
                p_degs[n-1] -= degree
            elif name.startswith("p"):
                n = int(name[1:])
                p_degs[n-1] += degree           
            else:
                # Unknown, we just pass it through (it should
                # correspond to a flavor group character).
                assert(name.startswith("F"))
                mon *= R.gen(i)^degree

    # The integration over the phases vanishes if there is any nonzero
    # dependence on the phases.
    if any(k != 0 for k in p_degs):
        c = 0

    return (c, mon)

# Take a weight vector, and return some canonical description valid as
# a Python identifier.
def canonical_irrep_var(weights):
    if all(w == 0 for w in weights):
        return 1
    return SR.var("Xi_"+"_".join([str(w) for w in weights]))

# Take a monomial in the characters and expand it as a sum of
# characters of irreps of the SU(3) flavor group. Quite similar to the
# integration over the gauge group above, really, only this time we
# keep all the information of the resulting representations. Returns a
# (linear) polynomial in the characters of the fundamental irreps of
# the flavor group.
def m_flavor_expand(m, LiE, setup_LiE=True):
    # Treat the constant monomial separately, since sage does not
    # treat it as an element of a polynomial ring.
    if m in QQ:
        return m

    R = m.parent()
    names = R.variable_names()

    # Set up LiE, we generically just need to do this once in
    # magnetic_integration, it is redundant to do so here for every
    # polynomial.
    if setup_LiE:
        LiE_AN_setup(2, LiE)

    # Split the monomial into terms, and extract the flavor rep.
    # This is a list of tuples (R,n,k). The total rep is the
    # tensor product of the individual reps. R denotes the irrep,
    # n is the degree of the Adams operator applied to it, and k
    # is the power to which to lift the resulting representation.
    rep = []
    for i, degree in enumerate(m.degrees()):
        if degree != 0:
            name = names[i]
            if name.startswith("Fc"):
                n = int(name[2:])
                rep.append(("afund", n, degree))
            elif name.startswith("F"):
                n = int(name[1:])
                rep.append(("fund", n, degree))
            else:
                # We should be able to recognize everything at
                # this point, error.
                print "Could not recognize '{0}'!".format(name)
                assert(1==0)

    Lrep = rep_to_LiE(rep)
    LiE.stdin.write(Lrep+"\n")
    result = 0
    for dim, irrep in LiE_read_module(LiE):
        if dim != 0:
            result += dim*canonical_irrep_var(irrep)

    return result

# Same as largeN_integration, but acts on the list of monomials
# produced by the plethystic exponential. It also does the flavor
# expansion.
def largeN_SOSU_and_flavor_integration(d, LiE, maxn,
                                       print_gauge_invariants=False):
    LiE_AN_setup(2, LiE)

    final = {}

    for p, mons in d:
        print "t^({0}) [{1} monomials] =>".format(p, len(mons)),
        if print_gauge_invariants:
            print
        L = []
        for J3, c, m in mons:
            k, mon = largeN_integration(m, maxn)
            if print_gauge_invariants:
                print "\t + ({3}) * {0} => {1}*{2}".format(m, k,
                                                             mon, c)
                print (J3, c*k, mon)

            if k != 0:
                L.append((J3, c*k, mon))

        if L == []:
            print "0"
            continue

        # Flavor expansion
        print "[{0}] =>".format(len(L)),
        partial = {}
        for J3, c, m in L:
            flavor_character = m_flavor_expand(m, LiE, setup_LiE=False)
            try:
                partial[(p,J3)] += c*flavor_character
            except KeyError:
                partial[(p,J3)] = c*flavor_character

        for charge, index in partial.iteritems():
            if index != 0:
                print "{0}: {1}".format(charge, index)
                final[charge] = index

        # In case all flavor indices cancelled
        if all(v == 0 for v in partial.itervalues()):
            print "0"

    return final

# Compute the large N SCI to the given degree in t (doing only the
# gauge integrations).
def largeN_I(degree, print_gauge_invariants=False):
    LiE = start_LiE()

    t = SR.var("t")

    maxn = int(degree*(3/2))
    names = []
    names += ["g{0}".format(n) for n in range(1,2*maxn+1)]
    names += ["r{0}".format(n) for n in range(1,2*maxn+1)]
    names += ["p{0}".format(n) for n in range(1,2*maxn+1)]
    names += ["pi{0}".format(n) for n in range(1,2*maxn+1)]
    names += ["F{0}".format(n) for n in range(1,maxn+1)]
    names += ["Fc{0}".format(n) for n in range(1,maxn+1)]

    # x is the chemical potential for J_3, and xi its inverse.
    names += ["x", "xi"]

    R = PolynomialRing(QQ, names=names)
    
    gs = R.gens()[:2*maxn]
    rs = R.gens()[2*maxn:4*maxn]
    ps = R.gens()[4*maxn:6*maxn]
    pis = R.gens()[6*maxn:8*maxn]
    Fs = R.gens()[8*maxn:9*maxn]
    Fcs = R.gens()[9*maxn:10*maxn]

    x, xi = R.gens()[-2:]

    t = SR.var("t")

    print "Building the plethystic exponent..."
    PE = sum(i_n(n, t^n, x^n, xi^n, gs, rs, ps, pis, Fs, Fcs)/n
             for n in range(1, maxn + 1))

    print "Starting formal Taylor expansion"
    T = expand(taylor(exp(PE), t, 0, degree))

    d = []
    print "Taylor expansion done, extracting terms"
    for coeff, p in T.coeffs(t):
        print "t^({0})".format(p)
        coeffR = R(coeff)
        L = []
        for c, m in zip(coeffR.coefficients(), coeffR.monomials()):
            degs = m.degrees()
            J3 = degs[-2] - degs[-1]
            L.append((J3, c, m(x=1,xi=1)))
        d.append((p, L))

    print "Integrating"
    d = largeN_SOSU_and_flavor_integration(d, LiE, maxn,
                                           print_gauge_invariants)
    print "Integration done"

    stop_LiE(LiE)

    return d

import subprocess as sp
import os.path
import re

# Open a pipe to LiE
def start_LiE(path=None):
    if path == None:
        path = "/home/inaki/Downloads/LiE"
    execname = os.path.join(path, "Lie.exe")
    try:
        LiE = sp.Popen([execname, "initfile", path],
                       close_fds=True,
                       stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
    except OSError as e:
        print "Got an error executing '{0}': {1}".format(execname,
                                                         e.strerror)
        return None
    return LiE

# Close a running LiE program
def stop_LiE(LiE):
    LiE.stdin.write("quit\n")
    LiE.wait()

# Representation of a single term.
def term_to_LiE(repname, n, k):
    if n > 1:
        repname = "Adams({0}, {1})".format(n, repname)
    if k > 1:
        repname = "p_tensor({0}, {1})".format(k, repname)
    return repname

# Take a representation in tensor product form, and construct the LiE
# representation (in the form of a string). This assumes that adjoint,
# fund and afund have been defined previously in the currently running
# LiE session.
def rep_to_LiE(rep):
    if len(rep) > 1:
        return "tensor({0}, {1})".format(term_to_LiE(*rep[0]),
                                         rep_to_LiE(rep[1:]))
    return term_to_LiE(*rep[0])

# Set up the LiE session for the A_N group.
def LiE_AN_setup(N, LiE, symm = False, asymm = False):
    LiE.stdin.write("setdefault A{0}\n".format(N))
    fundname = "1X"+str([1]+[0]*(N-1))
    LiE.stdin.write("fund = {0}\n".format(fundname))
    afundname = "1X"+str([0]*(N-1)+[1])
    LiE.stdin.write("afund = {0}\n".format(afundname))
    LiE.stdin.write("adj = adjoint()\n")
    if asymm:
        asymm = [0,1]+[0]*(N-2)
        LiE.stdin.write("asymm = 1X{0}\n".format(asymm))
        # complex conjugate of the antisymmetric
        asymmc = [0]*(N-2)+[1,0]
        LiE.stdin.write("asymmc = 1X{0}\n".format(asymmc))
    elif symm:
        symm = [2]+[0]*(N-1)
        LiE.stdin.write("symm = 1X{0}\n".format(symm))
        # complex conjugate of the symmetric
        symmc = [0]*(N-1)+[2]
        LiE.stdin.write("symmc = 1X{0}\n".format(symmc))

# Set up the LiE session for the A_N x A_N group.
def LiE_ANAN_setup(N, LiE):
    LiE.stdin.write("setdefault A{0}A{0}\n".format(N))
    fund1 = [1] + [0]*(2*N-1)
    LiE.stdin.write("fund1 = 1X{0}\n".format(fund1))
    afund1 = [0]*(N-1) + [1] + [0]*N
    LiE.stdin.write("afund1 = 1X{0}\n".format(afund1))
    fund2 = [0]*N + [1] + [0]*(N-1)
    LiE.stdin.write("fund2 = 1X{0}\n".format(fund2))
    afund2 = [0]*(2*N-1) + [1]
    LiE.stdin.write("afund2 = 1X{0}\n".format(afund2))
    bifund = [1] + [0]*(2*N-2) + [1] # Same as fund1 x afund2
    LiE.stdin.write("bifund = 1X{0}\n".format(bifund))
    abifund = [0]*(N-1) + [1,1] + [0]*(N-1)
    LiE.stdin.write("abifund = 1X{0}\n".format(abifund))

# Set up the LiE session for the B_N x A_M group. This assumes M > 1.
def LiE_BNAM_setup(N, M, LiE):
    LiE.stdin.write("setdefault B{0}A{1}\n".format(N,M))
    bifund = [1]+[0]*(N-1) + [0]*(M-1)+[1]
    LiE.stdin.write("bifund = 1X{0}\n".format(bifund))
    abifund = [1]+[0]*(N-1) + [1]+[0]*(M-1)
    LiE.stdin.write("abifund = 1X{0}\n".format(abifund))
    asymm = [0]*N + [0,1]+[0]*(M-2)
    LiE.stdin.write("asymm = 1X{0}\n".format(asymm))
    # complex conjugate of the antisymmetric
    asymmc = [0]*N + [0]*(M-2)+[1,0]
    LiE.stdin.write("asymmc = 1X{0}\n".format(asymmc))
    LiE.stdin.write("adj1 = adjoint()[2]\n")
    LiE.stdin.write("adj2 = adjoint()[1]\n")

# Set up the LiE session for the C_N x A_M group.
def LiE_CNAM_setup(N, M, LiE):
    LiE.stdin.write("setdefault C{0}A{1}\n".format(N,M))
    bifund = [1]+[0]*(N-1) + [0]*(M-1)+[1]
    LiE.stdin.write("bifund = 1X{0}\n".format(bifund))
    abifund = [1]+[0]*(N-1) + [1]+[0]*(M-1)
    LiE.stdin.write("abifund = 1X{0}\n".format(abifund))
    symm = [0]*N + [2]+[0]*(M-1)
    LiE.stdin.write("symm = 1X{0}\n".format(symm))
    # complex conjugate of the symmetric
    symmc = [0]*N + [0]*(M-1)+[2]
    LiE.stdin.write("symmc = 1X{0}\n".format(symmc))
    LiE.stdin.write("adj1 = adjoint()[2]\n")
    LiE.stdin.write("adj2 = adjoint()[1]\n")

# B1 = SO(3) is isomorphic to SU(2), and LiE treats it separately.
def _LiE_B1_setup(N, LiE):
    LiE.stdin.write("setdefault A1\n")
    LiE.stdin.write("fund = 1X[2]\n")
    LiE.stdin.write("adj = adjoint()\n")

def LiE_BN_setup(N, LiE):
    if N == 1:
        _LiE_B1_setup(N, LiE)
        return
    LiE.stdin.write("setdefault B{0}\n".format(N))
    fund = [1]+[0]*(N-1)
    LiE.stdin.write("fund = 1X{0}\n".format(fund))
    LiE.stdin.write("adj = adjoint()\n")

def LiE_CN_setup(N, LiE):
    LiE.stdin.write("setdefault C{0}\n".format(N))
    fund = [1]+[0]*(N-1)
    LiE.stdin.write("fund = 1X{0}\n".format(fund))
    LiE.stdin.write("adj = adjoint()\n")

def LiE_set_maxobjects(maxobjects, LiE):
    LiE.stdin.write("maxobjects {0}\n".format(maxobjects))

__LiE_pattern = re.compile(r"\s*([+-]?)\s*(\d+)X\[([\s\d,]*)]")

# Read a virtual module from LiE. The result is a list of (virtual
# dimension, irrep) tuples, where irrep is given as a list of integers.
def LiE_read_module(LiE):
    result = []
    while True:
        line = LiE.stdout.readline().strip()
        matches = re.findall(__LiE_pattern, line)
        if matches == []:
            print "Got '{0}', but didn't understand it, abort.".format(line)
            assert(1==0)
        for m in matches:
            dim = int(m[1])
            if m[0] == '-':
                dim = -dim
            weights = tuple(int(n) for n in m[2].replace(" ", "").split(","))
            result.append((dim, weights))
        if line[-1] != '+':
            break

    return result
