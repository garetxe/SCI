# Construct the letter for the electric theory given the chemical
# potentials. Here r denotes the R-charge of the quarks, and Xa and
# X{L,R} are the adjoint and matter characters respectively. F(a){L,R}
# are the flavor characters, possibly conjugated. Conventions are as
# in Dolan&Osborn.
def i_e(t, x, xi, q, qi, Xa, XL, XR, FL, FaL, FR, FaR, rq):
    P  = (2*t^2 - t*(x+xi))*Xa             # Vector multiplet
    P += t^rq*XL*FL*q - t^(2-rq)*XR*FaL*qi # left magnetic quark (q)
    P += t^rq*XR*FaR*qi - t^(2-rq)*XL*FR*q # right magnetic quark (q bar)
    P /= (1-t*x)*(1-t*xi)
    return P

# Electric plethystic exponential
def PE_e(bounds, t, Nf, Nc):
    rq = 1-Nc/Nf
    Bq = 1

    # Compute the maximum n that we need in the plethystic
    # exponential.
    maxn = 1
    while True:
        if all(deg > bounds[1] for deg in [rq*maxn, (2-rq)*maxn]):
            maxn -= 1
            break
        maxn += 1

    return _PE_em(bounds, t, Nf, Nc, rq, 0, Bq, maxn, meson = False)

# Similarly for the magnetic theory, with the difference that now we
# also have a dual meson. The conventions are also as in Dolan&Osborn.
def i_m(t, x, xi, q, qi, Xa, XL, XR,
        FL, FaL, FR, FaR, FM, FaM, rq, rM):
    P  = (2*t^2 - t*(x+xi))*Xa             # Vector multiplet
    P += t^rq*XL*FaL*q - t^(2-rq)*XR*FL*qi # left magnetic quark (q)
    P += t^rq*XR*FR*qi - t^(2-rq)*XL*FaR*q # right magnetic quark (q bar)
    P += t^rM*FM - t^(2-rM)*FaM            # dual meson
    P /= (1-t*x)*(1-t*xi)
    return P

# Magnetic plethystic exponential
def PE_m(bounds, t, Nf, Nc):
    rq = Nc/Nf
    rM = 2*(1-rq)
    Bq = Nc/(Nf-Nc)

    # Compute the maximum n that we need in the plethystic
    # exponential.
    maxn = 1
    while True:
        if all(deg > bounds[1] for deg in [rq*maxn, (2-rq)*maxn,
                                           rM*maxn, (2-rM)*maxn]):
            maxn -= 1
            break
        maxn += 1

    return _PE_em(bounds, t, Nf, Nc, rq, rM, Bq, maxn, meson = True)

# Much of the procedure is identical for the magnetic and electric
# computation, so this routine does the real job in both cases.
def _PE_em(bounds, t, Nf, Nc, rq, rM, Bq, maxn, meson):
    names = []
    names += ["XA%d"%(k+1,) for k in range(maxn)]
    names += ["XL%d"%(k+1,) for k in range(maxn)]
    names += ["XR%d"%(k+1,) for k in range(maxn)]
    names += ["FL%d"%(k+1,) for k in range(maxn)]
    names += ["FaL%d"%(k+1,) for k in range(maxn)]
    names += ["FR%d"%(k+1,) for k in range(maxn)] 
    names += ["FaR%d"%(k+1,) for k in range(maxn)]
    if meson:
        names += ["FM%d"%(k+1,) for k in range(maxn)]
        names += ["FaM%d"%(k+1,) for k in range(maxn)]

    # q is the baryonic chemical potential, and x the chemical
    # potential for J_3. qi and xi are their inverses.
    names += ["x", "xi", "q", "qi"]

    R = PolynomialRing(QQ, names=names)

    XAs = R.gens()[:maxn]
    XLs = R.gens()[maxn:2*maxn]
    XRs = R.gens()[2*maxn:3*maxn]
    FLs = R.gens()[3*maxn:4*maxn]
    FaLs = R.gens()[4*maxn:5*maxn]
    FRs = R.gens()[5*maxn:6*maxn]
    FaRs = R.gens()[6*maxn:7*maxn]
    if meson:
        FMs = R.gens()[7*maxn:8*maxn]
        FaMs = R.gens()[8*maxn:9*maxn]

    x, xi, q, qi = R.gens()[-4:]

    if meson:
        K = sum(i_m(t^k, x^k, xi^k, q^k, qi^k, XAs[k-1], XLs[k-1], XRs[k-1],
                    FLs[k-1], FaLs[k-1], FRs[k-1], FaRs[k-1],
                    FMs[k-1], FaMs[k-1], rq, rM)/k
                for k in range(1,maxn+1))
    else:
        K = sum(i_e(t^k, x^k, xi^k, q^k, qi^k, XAs[k-1], XLs[k-1], XRs[k-1],
                    FLs[k-1], FaLs[k-1], FRs[k-1], FaRs[k-1], rq)/k
                for k in range(1,maxn+1))

    # Taylor expand the result, coerce it into the polynomial ring,
    # and save it into a list for ease of manipulation later. Since
    # polynomial rings are by definition over positive powers of the
    # variables, we treat q and x separately (they can appear with
    # negative powers), and return a list of monomials,
    # with explicit B and J3 charge.
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
            B = (degs[-2] - degs[-1])*Bq
            J3 = degs[-4] - degs[-3]
            L.append((B, J3, c, m(x=1,xi=1,q=1,qi=1)))
        d.append((p, L))

    return d

# Given a monomial in the ring of characters above, returns the result
# of integrating it over the gauge group with the standard Haar
# measure (i.e., it counts the number of singlets in the associated
# rep under the gauge group).
def m_gauge_integration(m, N, LiE, setup_LiE=True):
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
        LiE_AN_setup(N-1, LiE)

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
            if name.startswith("XL"):
                n = int(name[2:])
                rep.append(("fund", n, degree))
            elif name.startswith("XR"):
                n = int(name[2:])
                rep.append(("afund", n, degree))
            elif name.startswith("XA"):
                n = int(name[2:])
                rep.append(("adj", n, degree))
            else:
                # Unknown, we just pass it through (it corresponds
                # to a flavor group character).
                assert(name.startswith("F"))
                mon *= R.gen(i)^degree

    return (n_singlets(N-1, rep, LiE), mon)

# Same as m_gauge_integration, but acts on the list of monomials
# produced by the plethystic exponential.
def gauge_integration(d, N, LiE):
    LiE_AN_setup(N-1, LiE)

    result = []
    for p, mons in d:
        L = []
        for B, J3, c, m in mons:
#            print ">>", m
            k,m = m_gauge_integration(m, N, LiE, setup_LiE=False)
#            print "<<", k, m
            if k != 0:
                L.append((B, J3, c*k, m))
        result.append((p, L))

    return result

# Take a weight vector, and return some canonical description valid as
# a Python identifier.
def canonical_irrep_var(weights):
    if all(w == 0 for w in weights):
        return 1
    return SR.var("Xi_"+"_".join([str(w) for w in weights]))

# Take a monomial in the characters and expand it as a sum of
# characters of irreps of the flavor group. Quite similar to the
# integration over the gauge group above, really, only this time we
# keep all the information of the resulting representations. Returns a
# (linear) polynomial in the characters of the fundamental irreps of
# the flavor group.
def m_flavor_expand(m, Nf, LiE, setup_LiE=True):
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
        LiE_ANAN_setup(Nf-1, LiE)

    # Split the monomial into terms, and extract the flavor rep.
    # This is a list of tuples (R,n,k). The total rep is the
    # tensor product of the individual reps. R denotes the irrep,
    # n is the degree of the Adams operator applied to it, and k
    # is the power to which to lift the resulting representation.
    rep = []
    for i, degree in enumerate(m.degrees()):
        if degree != 0:
            name = names[i]
            if name.startswith("FL"):
                n = int(name[2:])
                rep.append(("fund1", n, degree))
            elif name.startswith("FaL"):
                n = int(name[3:])
                rep.append(("afund1", n, degree))
            elif name.startswith("FR"):
                n = int(name[2:])
                rep.append(("fund2", n, degree))
            elif name.startswith("FaR"):
                n = int(name[3:])
                rep.append(("afund2", n, degree))
            elif name.startswith("FM"):
                n = int(name[2:])
                rep.append(("bifund", n, degree))
            elif name.startswith("FaM"):
                n = int(name[3:])
                rep.append(("abifund", n, degree))
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

# Take a dictionary, as produced by magnetic_integration, and expand
# the appearing characters in terms of irreps of the flavor group.
def flavor_expand(d2, Nf, LiE):
    LiE_ANAN_setup(Nf-1, LiE)

    result = {}
    for p, mons in d2:
        for B, J3, c, m in mons:
            flavor_character = m_flavor_expand(m, Nf, LiE, setup_LiE=False)
            try:
                result[(p,B,J3)] += c*flavor_character
            except KeyError:
                result[(p,B,J3)] = c*flavor_character

    # Clean up possible zero values, coming from cancellations in the
    # flavor character expansion.
    final = {}
    for charge, index in result.iteritems():
        if index != 0:
            final[charge] = index
    return final

# Compute the electric index in order n, with bounds[0] < n <= bounds[1] 
def electric_index(bounds, Nf, Nc):
    LiE = start_LiE()

    t = SR.var("t")
    d = PE_e(bounds, t, Nf, Nc)
    d2 = gauge_integration(d, Nc, LiE)
    print d2
    result = flavor_expand(d2, Nf, LiE)

    stop_LiE(LiE)

    return result

# Compute the magnetic index in order n, with bounds[0] < n <= bounds[1] 
def magnetic_index(bounds, Nf, Nc):
    LiE = start_LiE()

    t = SR.var("t")
    d = PE_m(bounds, t, Nf, Nc)
    d2 = gauge_integration(d, Nf-Nc, LiE)
    print d2
    result = flavor_expand(d2, Nf, LiE)

    stop_LiE(LiE)

    return result


# Returns how many singlets are there in the A_N representation
# given by rep.
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

# Flavor group representation (times 24) in the magnetic side at level 12/7:
# > p_tensor(4, fund1) + p_tensor(4, afund2) + 3*(p_tensor(2, Adams(2, fund1)) + p_tensor(2, Adams(2, afund2))) - 6*(tensor(p_tensor(2,fund1), Adams(2,fund1)) + tensor(p_tensor(2,afund2), Adams(2,afund2))) + 12*tensor(p_tensor(2, fund1), p_tensor(2, afund2)) + 8*( tensor(Adams(3, afund2), afund2) + tensor(Adams(3, fund1), fund1)) - 12*p_tensor(2, tensor(fund1, afund2)) + 12*tensor(Adams(2, fund1), Adams(2, afund2)) - 12*Adams(2, tensor(fund1, afund2)) - 6*Adams(4, fund1) - 6*Adams(4, afund2)
#     24X[0,0,0,0,0,0,0,0,1,0,0,0] +24X[0,0,0,1,0,0,0,0,0,0,0,0]

# (Vanishing) mesonic part only:
# 12*tensor(p_tensor(2, fund1), p_tensor(2, afund2)) - 12*p_tensor(2, tensor(fund1, afund2)) + 12*tensor(Adams(2, fund1), Adams(2, afund2)) - 12*Adams(2, tensor(fund1, afund2))

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
def LiE_AN_setup(N, LiE):
    LiE.stdin.write("setdefault A{0}\n".format(N))
    fundname = "1X"+str([1]+[0]*(N-1))
    LiE.stdin.write("fund = {0}\n".format(fundname))
    afundname = "1X"+str([0]*(N-1)+[1])
    LiE.stdin.write("afund = {0}\n".format(afundname))
    LiE.stdin.write("adj = adjoint()\n")

# Set up the LiE session for the A_NA_N group.
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
