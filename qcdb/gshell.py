"""Class to

"""


def INT_NCART(am):
    """Gives the number of cartesian functions for an angular momentum.
    define INT_NCART(am) ((am>=0) ? ((((am)+2)*((am)+1))>>1) : 0)

    """
    return (((am + 2) * (am + 1)) >> 1) if (am >= 0) else 0


def INT_NPURE(am):
    """Gives the number of spherical functions for an angular momentum.
    #define INT_NPURE(am) (2*(am)+1)

    """
    return 2 * am + 1


def INT_NFUNC(pu, am):
    """Gives the number of functions for an angular momentum based on pu.
    #define INT_NFUNC(pu,am) ((pu)?INT_NPURE(am):INT_NCART(am))

    """
    return INT_NPURE(am) if pu else INT_NCART(am)


def INT_CARTINDEX(am, i, j):
    """Computes offset index for cartesian function.
    #define INT_CARTINDEX(am,i,j) (((i) == (am))? 0 : (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)))

    """
    return 0 if (i == am) else ((((am - i + 1) * (am - i)) >> 1) + am - i - j)


def INT_ICART(a, b, c):
    """Given a, b, and c compute a cartesian offset.
    #define INT_ICART(a, b, c) (((((((a)+(b)+(c)+1)<<1)-(a))*((a)+1))>>1)-(b)-1)

    """
    return ((((((a + b + c + 1) << 1) - a) * (a + 1)) >> 1) - b - 1)


def INT_IPURE(l, m):
    """Given l and m compute a pure function offset.
    #define INT_IPURE(l, m) ((l)+(m))

    """
    return l + m


# Lookup array that when you index the angular momentum it returns the corresponding letter
PrimitiveType = ['Normalized', 'Unnormalized']
GaussianType = ['Cartesian', 'Pure']  # Cartesian = 0, Pure = 1


class ShellInfo(object):
    """This class has the same behavior as GaussianShell, but implements everything using
    slower data structures, which are easier to construct. These are used to build the
    basis set, which builds more efficient pointer-based GaussianShell objects.
    *  @param e An array of exponent values.
    *  @param am Angular momentum.
    *  @param pure Pure spherical harmonics, or Cartesian.
    *  @param c An array of contraction coefficients.
    *  @param nc The atomic center that this shell is located on. Must map
        back to the correct atom in the owning BasisSet molecule_. Used
        in integral derivatives for indexing.
    *  @param center The x, y, z position of the shell. This is passed to
        reduce the number of calls to the molecule.
    *  @param start The starting index of the first function this shell
        provides. Used to provide starting positions in matrices.
    *  @param pt Is the shell already normalized?

    """

    def __init__(self, am, c, e, pure, nc, center, start, pt='Normalized'):
        # Angular momentum
        self.l = am
        # Flag for pure angular momentum
        self.puream = pure
        # Exponents (of length nprimitives_)
        self.exp = e
        # Contraction coefficients (of length nprimitives_)
        self.coef = c
        # ERD normalized contraction coefficients (of length nprimitives_)
        self.erd_coef
        # Original (un-normalized) contraction coefficients (of length nprimitives)
        self.original_coef = [c[n] for n in range(len(c))]
        # Atom number this shell goes to. Needed when indexing integral derivatives.
        self.nc = nc
        # Atomic center number in the Molecule
        self.center = center
        #
        self.start = start
        # How many cartesian functions? (1=s, 3=p, 6=d, ...)
        self.ncartesian = INT_NCART(self.l)
        # How many functions? (1=s, 3=p, 5/6=d, ...) * Dependent on the value of puream_
        self.nfunction = INT_NFUNC(self.puream, self.l)

        # Compute the normalization constants
        if pt == 'Unnormalized':
            self.normalize_shell()
            self.erd_normalize_shell()

    def primitive_normalization(self, p):
        """Normalizes a single primitive.
        @param p The primitive index to normalize.
        @return Normalization constant to be applied to the primitive.

        """
        tmp1 = self.l + 1.5
        g = 2.0 * self.exp[p]
        z = pow(g, tmp1)
        return math.sqrt((pow(2.0, self.l) * z) / (math.pi * math.sqrt(math.pi) * df[2 * self.l]))

    def contraction_normalization(self):
        """Normalizes an entire contraction set. Applies the normalization to the coefficients
        *  @param gs The contraction set to normalize.

        """
        e_sum = 0.0
        for i in range(self.nprimitive()):
            for j in range(self.nprimitive()):
                g = self.exp[i] + self.exp[j]
                z = pow(g, self.l + 1.5)
                e_sum += self.coef[i] * self.coef[j] / z

        tmp = ((2.0 * math.pi / (2.0 / math.sqrt(math.pi))) * df[2 * self.l]) / pow(2.0, self.l)
        try:
            norm = math.sqrt(1.0 / (tmp * e_sum))
        except ZeroDivisionError:
            self.coef[i] = [1.0 for i in range(self.nprimitive())]
        # Set the normalization
        for i in range(self.nprimitive()):
            self.coef[i] *= norm

    def normalize_shell(self):
        """Handles calling primitive_normalization and
        contraction_normalization for you.

        """
        for i in range(self.nprimitive()):
            normalization = self.primitive_normalization(i)
            self.coef[i] *= normalization
        self.contraction_normalization()

    def erd_normalize_shell(self):
        """

        """
        self.erd_coef = []
        tsum = 0.0
        for j in range(self.nprimitive()):
            for k in range(j):
                a1 = self.exp[j]
                a2 = self.exp[k]
                temp = self.original_coef(j) * self.original_coef(k)
                temp2 = self.l + 1.5
                temp3 = 2.0 * math.sqrt(a1 * a2) / (a1 + a2)
                temp3 = pow(temp3, temp2)
                temp *= temp3
                tsum += temp
                if j != k:
                    tsum += temp
        prefac = 1.0
        if self.l > 1:
            prefac = pow(2.0, 2 * self.l) / df[2 * self.l]
        norm = math.sqrt(prefac / tsum)
        for j in range(self.nprimitive()):
            self.erd_coef.append(self.original_coef[j] * norm)

    def copy(self, nc=None, c=None):
        """Make a copy of the ShellInfo"""
        if nc is not None and c is not None:
            return ShellInfo(self.l, self.original_coef, self.exp,
                GaussianType[self.puream], nc, c,
                self.start, 'Unnormalized')
        else:
            return ShellInfo(self.l, self.original_coef, self.exp,
                GaussianType[self.puream], self.nc, self.center,
                self.start, 'Unnormalized')
        # better to just deepcopy?

    def nprimitive(self):
        """The number of primitive Gaussians"""
        return len(self.exp)

    def nfunction(self):
        """Total number of basis functions"""
        return INT_NFUNC(self.puream, self.l)

    def ncartesian(self):
        """Total number of functions if this shell was Cartesian"""
        return self.ncartesian

    def am(self):
        """The angular momentum of the given contraction"""
        return self.l

    def amchar(self):
        """The character symbol for the angular momentum of the given contraction"""
        return 'spdfghiklmnopqrtuvwxyz'[self.l]

    def AMCHAR(self):
        """The character symbol for the angular momentum of the given contraction (upper case)"""
        return 'SPDFGHIKLMNOPQRTUVWXYZ'[self.l]

    def is_cartesian(self):
        """Returns true if contraction is Cartesian"""
        return not self.puream

    def is_pure(self):
        """Returns true if contraction is pure"""
        return self.puream

    def center(self):
        """Returns the center of the Molecule this shell is on"""
        return self.center

    def ncenter(self):
        """Returns the atom number this shell is on. Used by integral derivatives for indexing."""
        return self.nc

    def exp(self, prim):
        """Returns the exponent of the given primitive"""
        return self.exp[prim]

    def coef(self, pi):
        """Return coefficient of pi'th primitive"""
        return self.coef[pi]

    def erd_coef(self, pi):
        """Return ERD normalized coefficient of pi'th primitive"""
        return self.erd_coef[pi]

    def original_coef(self, pi):
        """Return unnormalized coefficient of pi'th primitive"""
        return self.original_coef[pi]

    def exps(self):
        """Returns the exponent of the given primitive"""
        return self.exp

    def coefs(self):
        """Return coefficient of pi'th primitive and ci'th contraction"""
        return self.coef

    def original_coefs(self):
        """Return unnormalized coefficient of pi'th primitive and ci'th contraction"""
        return self.original_coef

    def print(self, outfile):
        """Print out the shell"""
        text = """    %c %3d 1.00\n""" % (AMCHAR(), self.nprimitive())
        for K in range(self.nprimitive()):
            text += """               %20.8f %20.8f\n""" % (self.exp[K], self.original_coef[K])
        with open(outfile, mode='w') as handle:
            handle.write(text)

    def normalize(self, l, m, n):
        """Normalize the angular momentum component"""
        return 1.0

    def function_index(self):
        """Basis function index where this shell starts."""
        return self.start

    def set_function_index(self, i):
        """Set basis function index where this shell starts."""
        self.start = i


#GaussianShell(0, nprimitive_,
#    uoriginal_coefficients_, ucoefficients_, uerd_coefficients_,
#    uexponents_, GaussianType(0), 0, xyz_, 0)
#
#GaussianShell(am, shell_nprim,
#    &uoriginal_coefficients_[ustart+atom_nprim], &ucoefficients_[ustart+atom_nprim], &uerd_coefficients_[ustart+atom_nprim],
#    &uexponents_[ustart+atom_nprim], puream, n, xyz_ptr, bf_count)
#
#GaussianShell(am, shell_nprim,
#    &uoriginal_coefficients_[prim_count], &ucoefficients_[prim_count], &uerd_coefficients_[prim_count],
#    &uexponents_[prim_count], puream, center, xyz_, bf_count)
#
#ShellInfo(am, contractions, exponents, gaussian_type, 0, center, 0, Unnormalized)
#
#    def constuctGS(self, am, nprimitive, oc, c, ec, e, pure, nc, center, start):
#        self.l = am
#        self.nprimitive = nprimitive
#        self.puream = pure
#        self.exp = e
#        self.original_coef = oc
#        self.coef = c
#        self.erd_coef = ec
#        self.nc = nc
#        self.center = center
#        self.start = start
#        self.ncartesian = INT_NCART(self.l)
#        self.nfunction = INT_NFUNC(self.puream, self.l)
