# sage: ISO = SOSU_index([2,3], N); ISO
# {(3, 1): -Xi_1_1 + Xi_2_2 - Xi_3_0 + Xi_6_0 + 3/2, (3, 3): 1/2, (3, -1): -Xi_1_1 + Xi_2_2 - Xi_3_0 + Xi_6_0 + 3/2, (3, -3): 1/2}
# sage: IUSp = USpSU_index([2,3], N-3); IUSp
# {(3, 1): -Xi_1_1 + Xi_2_2 - Xi_3_0 + Xi_6_0 + 1, (3, -1): -Xi_1_1 + Xi_2_2 - Xi_3_0 + Xi_6_0 + 1}

# Construct the letter for the G=(SO/USp)xSU theory given the chemical
# potentials. Here r denotes the R-charge of the quarks, and XA, Xa,
# Xb are characters under G of the adjoint and the reps associated
# with the fields A and B respectively. F is the character for the
# fundamental of the SU(3) flavor group. A "c" at the end denotes
# complex conjugation. We also introduce xi = 1/x. Other conventions
# are as in Dolan&Osborn.
def i_SOUSp(maxdegree, k, t, x, xi, XA, Xa, Xac, Xb, Xbc, F, Fc, ra, rb):
    P = 0
    # P  = (2*t^2 - t*(x+xi))*XA             # Vector multiplet
    if 2*k <= maxdegree:
        P += 2*t^(2*k)*XA
    if k <= maxdegree:
        P -= t^k*(x+xi)*XA
    # P += t^ra*Xa*F - t^(2-ra)*Xac*Fc       # A (the bifundamental)
    if k*ra <= maxdegree:
        P += t^(k*ra)*Xa*F
    if k*(2-ra) <= maxdegree:
        P -= t^(k*(2-ra))*Xac*Fc
    # P += t^rb*Xb*F - t^(2-rb)*Xbc*Fc       # B ((anti)symmetric)
    if k*rb <= maxdegree:
        P += t^(k*rb)*Xb*F
    if k*(2-rb) <= maxdegree:
        P -= t^(k*(2-rb))*Xbc*Fc

    P /= (1-t^k*x)*(1-t^k*xi)

    return taylor(P, t, 0, ceil(maxdegree))

# An older version, obviously correct.
def _i_SOUSp(t, x, xi, XA, Xa, Xac, Xb, Xbc, F, Fc, ra, rb):
    P  = (2*t^2 - t*(x+xi))*XA             # Vector multiplet
    P += t^ra*Xa*F - t^(2-ra)*Xac*Fc       # A (the bifundamental)
    P += t^rb*Xb*F - t^(2-rb)*Xbc*Fc       # B ((anti)symmetric)
    P /= (1-t*x)*(1-t*xi)
    return P


# Plethystic exponential for the SO(N-4)xSU(N) group.
def PE_SOSU(bounds, N):
    ra = 2/3 + 2/N
    rb = 2/3 - 4/N

    return _PE_SOUSp(bounds, ra, rb)

# Plethystic exponential for the USp(N+4)xSU(N) group.
def PE_USpSU(bounds, N):
    ra = 2/3 - 2/N
    rb = 2/3 + 4/N

    return _PE_SOUSp(bounds, ra, rb)

# The plethystic exponential itself is identical between the SO and
# USp cases.
def _PE_SOUSp(bounds, ra, rb):
    # Compute the maximum n that we need in the plethystic
    # exponential to expand up to the given order.
    maxn = 1
    while True:
        if all(deg > bounds[1] for deg in [maxn,
                                           ra*maxn, (2-ra)*maxn,
                                           rb*maxn, (2-rb)*maxn]):
            maxn -= 1
            break
        maxn += 1

    print "maxn is", maxn
    print "ra", ra
    print "rb", rb

    names = []
    # Adjoint of the B/C factor of the gauge group
    names += ["XABC%d"%(k+1,) for k in range(maxn)]
    # Adjoint of the A factor
    names += ["XAA%d"%(k+1,) for k in range(maxn)]
    names += ["Xa%d"%(k+1,) for k in range(maxn)]
    names += ["Xac%d"%(k+1,) for k in range(maxn)]
    names += ["Xb%d"%(k+1,) for k in range(maxn)]
    names += ["Xbc%d"%(k+1,) for k in range(maxn)]
    names += ["F%d"%(k+1,) for k in range(maxn)]
    names += ["Fc%d"%(k+1,) for k in range(maxn)]

    # x is the chemical potential for J_3, and xi its inverse.
    names += ["x", "xi"]

    R = PolynomialRing(QQ, names=names)

    XABCs = R.gens()[:maxn]
    XAAs = R.gens()[maxn:2*maxn]
    Xas = R.gens()[2*maxn:3*maxn]
    Xacs = R.gens()[3*maxn:4*maxn]
    Xbs = R.gens()[4*maxn:5*maxn]
    Xbcs = R.gens()[5*maxn:6*maxn]
    Fs = R.gens()[6*maxn:7*maxn]
    Fcs = R.gens()[7*maxn:8*maxn]

    x, xi = R.gens()[-2:]

    t = SR.var("t")

    print "Building the plethystic exponent..."
    K = 0
    for k in range(1, maxn+1):
        K += i_SOUSp(bounds[1], k, t, x^k, xi^k, (XABCs[k-1] + XAAs[k-1]),
                     Xas[k-1], Xacs[k-1], Xbs[k-1], Xbcs[k-1],
                     Fs[k-1], Fcs[k-1], ra, rb)/k

    print K

    d = []
    print "Starting formal Taylor expansion"
    T = taylor(exp(K), t, 0, bounds[1])

    print "Taylor expansion done, extracting terms"
    for coeff, p in T.coeffs(t):
        if p <= bounds[0]:
            continue
        print "t^({0})".format(p)
        coeffR = R(coeff)
        L = []
        for c, m in zip(coeffR.coefficients(), coeffR.monomials()):
            degs = m.degrees()
            J3 = degs[-2] - degs[-1]
            L.append((J3, c, m(x=1,xi=1)))
        d.append((p, L))

    print "Done with the plethystic exponential."

    return d

# And older version, more obviously correct.
def __PE_SOUSp(bounds, t, ra, rb):
    # Compute the maximum n that we need in the plethystic
    # exponential to expand up to the given order.
    maxn = 1
    while True:
        if all(deg > bounds[1] for deg in [maxn,
                                           ra*maxn, (2-ra)*maxn,
                                           rb*maxn, (2-rb)*maxn]):
            maxn -= 1
            break
        maxn += 1
    print "maxn is", maxn
    print "ra", ra
    print "rb", rb

    names = []
    # Adjoint of the B/C factor of the gauge group
    names += ["XABC%d"%(k+1,) for k in range(maxn)]
    # Adjoint of the A factor
    names += ["XAA%d"%(k+1,) for k in range(maxn)]
    names += ["Xa%d"%(k+1,) for k in range(maxn)]
    names += ["Xac%d"%(k+1,) for k in range(maxn)]
    names += ["Xb%d"%(k+1,) for k in range(maxn)]
    names += ["Xbc%d"%(k+1,) for k in range(maxn)]
    names += ["F%d"%(k+1,) for k in range(maxn)]
    names += ["Fc%d"%(k+1,) for k in range(maxn)]

    # x is the chemical potential for J_3, and xi its inverse.
    names += ["x", "xi"]

    R = PolynomialRing(QQ, names=names)

    XABCs = R.gens()[:maxn]
    XAAs = R.gens()[maxn:2*maxn]
    Xas = R.gens()[2*maxn:3*maxn]
    Xacs = R.gens()[3*maxn:4*maxn]
    Xbs = R.gens()[4*maxn:5*maxn]
    Xbcs = R.gens()[5*maxn:6*maxn]
    Fs = R.gens()[6*maxn:7*maxn]
    Fcs = R.gens()[7*maxn:8*maxn]

    x, xi = R.gens()[-2:]

    K = sum(i_SOUSp(t^k, x^k, xi^k, (XABCs[k-1] + XAAs[k-1]),
                    Xas[k-1], Xacs[k-1], Xbs[k-1], Xbcs[k-1],
                    Fs[k-1], Fcs[k-1], ra, rb)/k
            for k in range(1,maxn+1))

    # Taylor expand the result, coerce it into the polynomial ring,
    # and save it into a list for ease of manipulation later. Since
    # polynomial rings are by definition over positive powers of the
    # variables, we treat q and x separately (they can appear with
    # negative powers), and return a list of monomials,
    # with explicit J3 charge.
    d = []
    print "Starting formal Taylor expansion"
    T = taylor(exp(K), t, 0, ceil(bounds[1]))
    print "Taylor expansion done, extracting terms"
    for coeff, p in T.coeffs(t):
        if p <= bounds[0] or p > bounds[1]:
            continue
        print "t^({0})".format(p)
        coeffR = R(coeff)
        L = []
        for c, m in zip(coeffR.coefficients(), coeffR.monomials()):
            degs = m.degrees()
            J3 = degs[-2] - degs[-1]
            L.append((J3, c, m(x=1,xi=1)))
        d.append((p, L))

    return d

# Given a monomial in the ring of characters above, returns the result
# of integrating it over the SO(N-4)xSU(N) gauge group with the
# standard Haar measure (i.e., it counts the number of singlets in the
# associated rep under the gauge group). Here N is expected to be
# odd. We integrate first over SU(N) and then over SO(N-4).
def m_SOSU_integration(m, N, LiE):
    # Treat the constant monomial separately, since sage does not
    # treat it as an element of a polynomial ring.
    if m in QQ:
        return (1, m)

    R = m.parent()
    names = R.variable_names()

    LiE_AN_setup(N-1, LiE, asymm=True)

    mon = 1

    # Split the monomial into terms, and extract the gauge rep.
    # This is a list of tuples (R,n,k). The total rep is the
    # tensor product of the individual reps. R denotes the irrep,
    # n is the degree of the Adams operator applied to it, and k
    # is the power to which to lift the resulting representation.
    rep = []
    for i, degree in enumerate(m.degrees()):
        if degree != 0:
            name = names[i]
            if name.startswith("XAA"):
                n = int(name[3:])
                rep.append(("adj", n, degree))
            elif name.startswith("Xac"):
                n = int(name[3:])
                rep.append(("fund", n, degree))
                # This is the bifundamental, keep it for the SO(N)
                # integration.
                mon *= R.gen(i)^degree
            elif name.startswith("Xa"):
                n = int(name[2:])
                rep.append(("afund", n, degree))
                mon *= R.gen(i)^degree
            elif name.startswith("Xbc"):
                n = int(name[3:])
                rep.append(("asymmc", n, degree))
            elif name.startswith("Xb"):
                n = int(name[2:])
                rep.append(("asymm", n, degree))
            else:
                # Unknown, we just pass it through (it corresponds
                # to a flavor or SO group character).
                mon *= R.gen(i)^degree

    k = n_singlets(N-1, rep, LiE)
    if k == 0:
        return (0, 1)

    # Got something non-vanishing, integrate over the SO factor
    m = mon
    mon = 1

    if m in QQ:
        return (k, m)
    
    r = (N-5)/2
    assert (r in ZZ)
    LiE_BN_setup(int(r), LiE)

    rep = []
    for i, degree in enumerate(m.degrees()):
        if degree != 0:
            name = names[i]
            if name.startswith("XABC"):
                n = int(name[4:])
                rep.append(("adj", n, degree))
            elif name.startswith("Xac"):
                n = int(name[3:])
                rep.append(("fund", n, degree))
            elif name.startswith("Xa"):
                n = int(name[2:])
                rep.append(("fund", n, degree))
            else:
                # Unknown, we just pass it through (it must correspond
                # to a flavor group character).
                assert(name.startswith("F"))
                mon *= R.gen(i)^degree

    return (k*n_singlets(int(r), rep, LiE), mon)

# Older, slower but more obviously correct version of the above, for
# comparison.
def _m_SOSU_integration(m, N, LiE, setup_LiE=True):
    # Treat the constant monomial separately, since sage does not
    # treat it as an element of a polynomial ring.
    if m in QQ:
        return (1, m)

    R = m.parent()
    names = R.variable_names()

    # Set up LiE, we generically just need to do this once in
    # magnetic_integration, it is redundant to do so here for every
    # polynomial.
    if setup_LiE:
        r = (N-5)/2
        assert (r in ZZ)
        LiE_BNAM_setup(int(r), N-1, LiE)

    # The remaining monomial, after integrating over the gauge
    # group. It will be a funcion of the flavor characters only.
    mon = 1

    # Split the monomial into terms, and extract the gauge rep.
    # This is a list of tuples (R,n,k). The total rep is the
    # tensor product of the individual reps. R denotes the irrep,
    # n is the degree of the Adams operator applied to it, and k
    # is the power to which to lift the resulting representation.
    rep = []
    for i, degree in enumerate(m.degrees()):
        if degree != 0:
            name = names[i]
            if name.startswith("XABC"):
                n = int(name[4:])
                rep.append(("adj1", n, degree))
            elif name.startswith("XAA"):
                n = int(name[3:])
                rep.append(("adj2", n, degree))
            elif name.startswith("Xac"):
                n = int(name[3:])
                rep.append(("abifund", n, degree))
            elif name.startswith("Xa"):
                n = int(name[2:])
                rep.append(("bifund", n, degree))
            elif name.startswith("Xbc"):
                n = int(name[3:])
                rep.append(("asymmc", n, degree))
            elif name.startswith("Xb"):
                n = int(name[2:])
                rep.append(("asymm", n, degree))
            else:
                # Unknown, we just pass it through (it corresponds
                # to a flavor group character).
                assert(name.startswith("F"))
                mon *= R.gen(i)^degree

    return (n_singlets(int((3*N-7)/2), rep, LiE), mon)

# Same as m_SOSU_integration, but acts on the list of monomials
# produced by the plethystic exponential. It also does the flavor
# expansion.
def SOSU_and_flavor_integration(d, N, LiE, print_gauge_invariants=False):
    final = {}
    for p, mons in d:
        print "t^({0}) [{1} monomials] =>".format(p, len(mons)),
        if print_gauge_invariants:
            print
        L = []
        for J3, c, m in mons:
            k,mon = m_SOSU_integration(m, N, LiE)
            if print_gauge_invariants and J3 == 2:
                print "\t + ({3}) * {0} => {1}*{2}".format(m, k,
                                                             mon, c)
                print (J3, c*k, mon)

            if k != 0:
                L.append((J3, c*k, mon))
        if L == []:
            print "0"
            continue

        # Flavor expansion
        LiE_AN_setup(2, LiE)
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

# Given a monomial in the ring of characters above, returns the result
# of integrating it over the USp(N+4)xSU(N) gauge group with the
# standard Haar measure (i.e., it counts the number of singlets in the
# associated rep under the gauge group). Here N is expected to be
# even.
def m_USpSU_integration(m, N, LiE, setup_LiE=True):
    # Treat the constant monomial separately, since sage does not
    # treat it as an element of a polynomial ring.
    if m in QQ:
        return (1, m)

    R = m.parent()
    names = R.variable_names()

    # Set up LiE, we generically just need to do this once in
    # magnetic_integration, it is redundant to do so here for every
    # polynomial.
    if setup_LiE:
        r = (N+4)/2
        assert (r in ZZ)
        LiE_CNAM_setup(int(r), N-1, LiE)

    # The remaining monomial, after integrating over the gauge
    # group. It will be a funcion of the flavor characters only.
    mon = 1

    # Split the monomial into terms, and extract the gauge rep.
    # This is a list of tuples (R,n,k). The total rep is the
    # tensor product of the individual reps. R denotes the irrep,
    # n is the degree of the Adams operator applied to it, and k
    # is the power to which to lift the resulting representation.
    rep = []
    for i, degree in enumerate(m.degrees()):
        if degree != 0:
            name = names[i]
            if name.startswith("XAA"):
                n = int(name[3:])
                rep.append(("adj2", n, degree))
            elif name.startswith("XABC"):
                n = int(name[4:])
                rep.append(("adj1", n, degree))
            elif name.startswith("Xac"):
                n = int(name[3:])
                rep.append(("abifund", n, degree))
            elif name.startswith("Xa"):
                n = int(name[2:])
                rep.append(("bifund", n, degree))
            elif name.startswith("Xbc"):
                n = int(name[3:])
                rep.append(("symmc", n, degree))
            elif name.startswith("Xb"):
                n = int(name[2:])
                rep.append(("symm", n, degree))
            else:
                # Unknown, we just pass it through (it corresponds
                # to a flavor group character).
                assert(name.startswith("F"))
                mon *= R.gen(i)^degree

    return (n_singlets(int((3*N+2)/2), rep, LiE), mon)

# Same as m_SOSU_integration, but for the dual theory.
def USpSU_and_flavor_integration(d, N, LiE, print_gauge_invariants=False):
    r = (N+4)/2
    assert (r in ZZ), "USp rank not an integer!"

    final = {}
    for p, mons in d:
        print "t^({0}) [{1} monomials] =>".format(p, len(mons)),
        if print_gauge_invariants:
            print
        LiE_CNAM_setup(int(r), N-1, LiE)
        L = []
        for J3, c, m in mons:
            k,mon = m_USpSU_integration(m, N, LiE, setup_LiE=False)
            if k != 0:
                if print_gauge_invariants:
                    print "\t + {0} => {1}*{2}".format(m, k, mon)
                L.append((J3, c*k, mon))

        if L == []:
            print "0"
            continue

        # Flavor integration
        LiE_AN_setup(2, LiE)
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

# Compute the SO(N-4)xSU(N) index to order n, with bounds[0] < n <=
# bounds[1]
def SOSU_index(bounds, N, print_gauge_invariants=False):
    LiE = start_LiE()

    LiE_set_maxobjects(19999999, LiE)

    d = PE_SOSU(bounds, N)
    print "Gauge integration and flavor expansion..."
    result = SOSU_and_flavor_integration(d, N, LiE, print_gauge_invariants)

    stop_LiE(LiE)

    return result

# Compute the USp(N+4)xSU(N) index to order n, with bounds[0] < n <=
# bounds[1]
def USpSU_index(bounds, N, print_gauge_invariants=False):
    LiE = start_LiE()

    LiE_set_maxobjects(3999999, LiE)

    d = PE_USpSU(bounds, N)
    print "Gauge integration and flavor expansion..."
    result = USpSU_and_flavor_integration(d, N, LiE, print_gauge_invariants)

    stop_LiE(LiE)

    return result

# Returns how many singlets are there in the given representation. N
# is the rank of the group.
def n_singlets(N, rep, LiE):
    # The trivial representation
    if rep == []:
        return 1

    Lrep = rep_to_LiE(rep)
    LiE.stdin.write(Lrep+"\n")

    mod = LiE_read_module(LiE)
    trivial = tuple([0]*N)
    trivial_reps = [k for k, R in mod if R == trivial]
    assert (len(trivial_reps) < 2)

    if len(trivial_reps) == 0:
        return 0

    return trivial_reps[0]


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
