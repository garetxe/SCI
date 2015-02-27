# Given a partition with a maximum of 3 rows (corresponding to a rep
# of SU(3)), return the Dynkin labels for the representation, in LiE's
# notation. This is essentially the same as doing plethysm(partition,
# 1X[1,0]) in LiE.
#
# For instance [10,3,2] -> [7,1]
def partition_to_dynkin(p):
    if (len(p) > 3):
        raise NotImplementedError, "Partitions with more than 3 rows not supported"
    # Normalize the partition to exactly two rows.
    if len(p) == 3:
        p = [p[0]-p[2], p[1]-p[2]]
    elif len(p) == 1:
        p = [p[0], 0]
    elif len(p) == 0:
        p = [0,0]

    return [p[0]-p[1], p[1]]

def B21():
    # Partitions with more than 3 rows vanish on the flavor side.
    partitions = [part for part in partitions_list(21) if len(part) < 4]
    print "There are {0} partitions to check".format(len(partitions))

    N = 6

    nvparts = []

    for part in partitions:
        # Here "." refers to the path where the LiE executable is located,
        # in this case the code expects it to be in the current running
        # directory.
        LiE = start_LiE(".")

        LiE_set_maxobjects(39999999, LiE)
        LiE_AN_setup(N, LiE, asymm = True)
        
        n = n_singlets(N, LiE_plethysm(part, LiE))
        print "{0}{1}=> {2}".format(part, " "*(12-len(str(part))), n)

        if n != 0:
            nvparts.append((n, part))

        stop_LiE(LiE)

    # Write down the result in LiE notation.
    rep = " + ".join("{0}X{1}".format(c, partition_to_dynkin(part))
                     for c,part in nvparts)

    print "Got the following result:"
    print rep

    return rep

import subprocess as sp
import os.path
import re

# Return the module for the given plethysm, obtained with LiE
def LiE_plethysm(partition, LiE):
    LiE.stdin.write("plethysm({0}, asymm)\n".format(partition))

    return LiE_read_module(LiE)

# Returns how many singlets are there in the given module read from
# LiE, for a group with rank N.
def n_singlets(N, mod):
    trivial = tuple([0]*N)
    trivial_reps = [k for k, R in mod if R == trivial]
    assert (len(trivial_reps) < 2)

    if len(trivial_reps) == 0:
        return 0

    return trivial_reps[0]

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
        # This is a debug message, we can ignore
        if line.startswith("[CD]"):
            if not line[4:].strip().startswith("Computing"):
                # Things being added to the cache is interesting
                print line
            continue

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
