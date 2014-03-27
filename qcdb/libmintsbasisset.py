import os
import re
import collections
from exceptions import *
from psiutil import search_file
from molecule import Molecule
from libmintsgshell import GaussianShell
from libmintsbasissetparser import Gaussian94BasisSetParser


class BasisSet(object):
    """Basis set container class
    Reads the basis set from a checkpoint file object. Also reads the molecule
    from the checkpoint file storing the information in an internal Molecule class
    which can be accessed using molecule().

    """

    # <<< Globals >>>

    # Has static information been initialized?
    initialized_shared = False
    # Global arrays of x, y, z exponents (Need libmint for max ang mom)
    LIBINT_MAX_AM = 6  # TODO
    exp_ao = [[] for l in range(LIBINT_MAX_AM)]

    def __init__(self, *args):

        # <<< Basis BasisSet Information >>>

        # The name of this basis set (e.g. "BASIS", "RI BASIS")
        self.name = None
        # Array of gaussian shells
        self.shells = None
        # Molecule object.
        self.molecule = None

        # <<< Scalars >>>

        # Number of atomic orbitals (Cartesian)
        self.PYnao = None
        # Number of basis functions (either cartesian or spherical)
        self.PYnbf = None
        # The number of unique primitives
        self.n_uprimitive = None
        # The number of shells
        self.n_shells = None
        # The number of primitives
        self.PYnprimitive = None
        # The maximum angular momentum
        self.PYmax_am = None
        # The maximum number of primitives in a shell
        self.max_nprimitive = None
        # Whether the basis set is uses spherical basis functions or not
        self.puream = None

        # <<< Arrays >>>

        # The number of primitives (and exponents) in each shell
        self.n_prim_per_shell = None
        # The first (Cartesian) atomic orbital in each shell
        self.shell_first_ao = None
        # The first (Cartesian / spherical) basis function in each shell
        self.shell_first_basis_function = None
        # Shell number to atomic center.
        self.shell_center = None
        # Which shell does a given (Cartesian / spherical) function belong to?
        self.function_to_shell = None
        # Which shell does a given Cartesian function belong to?
        self.ao_to_shell = None
        # Which center is a given function on?
        self.function_center = None
        # How many shells are there on each center?
        self.center_to_nshell = None
        # What's the first shell on each center?
        self.center_to_shell = None

        # The flattened lists of unique exponents
        self.uexponents = None
        # The flattened lists of unique contraction coefficients (normalized)
        self.ucoefficients = None
        # The flattened lists of unique contraction coefficients (as provided by the user)
        self.uoriginal_coefficients = None
        # The flattened lists of ERD normalized contraction coefficients
        self.uerd_coefficients = None
        # The flattened list of Cartesian coordinates for each atom
        self.xyz = None

        for item in args:
            print type(item)
        if len(args) == 0:
            self.constructor_zero_ao_basis()
        elif len(args) == 3 and \
            isinstance(args[0], basestring) and \
            isinstance(args[1], Molecule) and \
            isinstance(args[2], collections.OrderedDict):
            #isinstance(args[0], Gaussian94BasisSetParser) and \
            #isinstance(args[1], Molecule) and \
            #isinstance(args[2], basestring):
            self.constructor_role_mol_shellmap(*args)
        elif len(args) == 2:
            # TODO typecheck
            self.constructor_basisset_center(args)
        else:
            raise ValidationError('BasisSet::constructor: Inappropriate configuration of constructor arguments')

    def constructor_zero_ao_basis(self):
        """Constructs a zero AO basis set"""

        if not self.initialized_shared:
            self.initialize_singletons()
        self.initialized_shared = True

        # Add a dummy atom at the origin, to hold this basis function
        self.molecule = Molecule()
        self.molecule.add_atom(0, 0.0, 0.0, 0.0)
        # Fill with data representing a single S function, at the origin, with 0 exponent
        self.n_uprimitive = 1
        self.n_shells = 1
        self.PYnprimitive = 1
        self.PYnao = 1
        self.PYnbf = 1
        self.uerd_coefficients = [1.0]
        self.n_prim_per_shell = [1]
        self.uexponents = [0.0]
        self.ucoefficients = [1.0]
        self.uoriginal_coefficients = [1.0]
        self.shell_first_ao = [0]
        self.shell_first_basis_function = [0]
        self.ao_to_shell = [0]
        self.function_to_shell = [0]
        self.function_center = [0]
        self.shell_center = [0]
        self.center_to_nshell = [0]
        self.center_to_shell = [0]
        self.puream = False
        self.PYmax_am = 0
        self.max_nprimitive = 1
        self.xyz = [0.0, 0.0, 0.0]
        self.name = '(Empty Basis Set)'
        self.shells = []
        self.shells.append(GaussianShell(0, self.PYnprimitive,
            self.uoriginal_coefficients, self.ucoefficients, self.uerd_coefficients,
            self.uexponents, 'Cartesian', 0, self.xyz, 0))

    def constructor_role_mol_shellmap(self, role, mol, shell_map):
        """

        """
        self.molecule = mol
        self.name = role

        # Singletons
        if not self.initialized_shared:
            self.initialize_singletons()
        self.initialized_shared = True

        natom = self.molecule.natom()

        # These will tell us where the primitives for [basis][symbol] start and end, in the compact array
        #std::map<std::string, std::map<std::string, int > >  primitive_start;
        #std::map<std::string, std::map<std::string, int > >  primitive_end;
        primitive_start = {}
        primitive_end = {}

        # First, loop over the unique primitives, and store them
        uexps = []
        ucoefs = []
        uoriginal_coefs = []
        uerd_coefs = []
        self.n_uprimitive = 0
        for basisfirst, basissecond in shell_map.items():
            basis = basisfirst
            symbol_map = shell_map[basis]
            primitive_start[basis] = {}
            primitive_end[basis] = {}
            for symbolfirst, symbolsecond in symbol_map.items():
                #symbol = symbolfirst
                #shells = symbol_map[symbol]
                #primitive_start[basis][symbol] = self.n_uprimitive
                label = symbolfirst
                shells = symbol_map[label]
                primitive_start[basis][label] = self.n_uprimitive
                for i in range(len(shells)):
                    shell = shells[i]
                    #print shell
                    for prim in range(shell.nprimitive()):
                        uexps.append(shell.exp(prim))
                        #print basisfirst, symbolfirst, i, prim, shell.exp(prim), len(uexps)
                        ucoefs.append(shell.coef(prim))
                        uoriginal_coefs.append(shell.original_coef(prim))
                        uerd_coefs.append(shell.erd_coef(prim))
                        self.n_uprimitive += 1
                #primitive_end[basis][symbol] = self.n_uprimitive
                primitive_end[basis][label] = self.n_uprimitive

        #print 'prim_stt', primitive_start
        #print 'prim_end', primitive_end

        # Count basis functions, shells and primitives
        self.n_uprimitive = 0
        self.n_shells = 0
        self.PYnprimitive = 0
        self.PYnao = 0
        self.PYnbf = 0
        for n in range(natom):
            atom = self.molecule.atom_entry(n)
            basis = atom.basisset(role)
            #symbol = atom.symbol()
            #shells = shell_map[basis][symbol]
            label = atom.label()
            shells = shell_map[basis][label]
            for i in range(len(shells)):
                shell = shells[i]
                nprim = shell.nprimitive()
                self.n_uprimitive += nprim
                self.PYnprimitive += nprim
                self.n_shells += 1
                self.PYnao += shell.ncartesian()
                self.PYnbf += shell.nfunction()
                #print 'PYnprim', self.PYnprimitive, 'n', n, 'i', i, 'nprim', nprim
                #print "natom %d, n_uprim %d, nprim %d, nao %d, nbf %d, n_shells %d\n" % \
                #    (natom, self.n_uprimitive, self.PYnprimitive, self.PYnao, self.PYnbf, self.n_shells)

        # Allocate arrays
        self.n_prim_per_shell = [0] * self.n_shells
        # The unique primitives
        self.uexponents = [0.0] * self.n_uprimitive
        self.ucoefficients = [0.0] * self.n_uprimitive
        self.uoriginal_coefficients = [0.0] * self.n_uprimitive
        self.uerd_coefficients = [0.0] * self.n_uprimitive
        for i in range(len(uexps)):  # TODO change from libmints
        #for i in range(self.n_uprimitive):
            #print i, self.n_uprimitive, len(self.uexponents), len(uexps)
            self.uexponents[i] = uexps[i]
            self.ucoefficients[i] = ucoefs[i]
            self.uoriginal_coefficients[i] = uoriginal_coefs[i]
            self.uerd_coefficients[i] = uerd_coefs[i]

        self.shell_first_ao = [0] * self.n_shells
        self.shell_first_basis_function = [0] * self.n_shells
        self.shells = [None] * self.n_shells
        self.ao_to_shell = [0] * self.PYnao
        self.function_to_shell = [0] * self.PYnbf
        self.function_center = [0] * self.PYnbf
        self.shell_center = [0] * self.n_shells
        self.center_to_nshell = [0] * natom
        self.center_to_shell = [0] * natom
        self.xyz = [0.0] * 3 * natom

        # Now loop over all atoms, and point to the appropriate unique data
        shell_count = 0
        ao_count = 0
        bf_count = 0
        xyz_ptr = self.xyz  # TODO
        self.puream = False
        self.PYmax_am = 0
        self.max_nprimitive = 0
        for n in range(natom):
            atom = self.molecule.atom_entry(n)
            basis = atom.basisset(role)
            #symbol = atom.symbol()
            #shells = shell_map[basis][symbol]
            #ustart = primitive_start[basis][symbol]
            #uend = primitive_end[basis][symbol]
            label = atom.label()
            shells = shell_map[basis][label]
            ustart = primitive_start[basis][label]
            uend = primitive_end[basis][label]
            nshells = len(shells)
            self.center_to_nshell[n] = nshells
            self.center_to_shell[n] = shell_count
            atom_nprim = 0
            for i in range(nshells):
                thisshell = shells[i]
                self.shell_first_ao[shell_count] = ao_count
                self.shell_first_basis_function[shell_count] = bf_count
                shell_nprim = thisshell.nprimitive()
                am = thisshell.am()
                self.max_nprimitive = max(shell_nprim, self.max_nprimitive)
                self.PYmax_am = max(am, self.PYmax_am)
                self.shell_center[shell_count] = n
                puream = 'Pure' if thisshell.is_pure() else 'Cartesian'
                self.puream = thisshell.is_pure()
                #if puream:
                #    self.puream = True
                self.puream = thisshell.is_pure()
                #print 'PUREAM: ', puream, self.puream
                #print "atom %d basis %s shell %d nprim %d atom_nprim %d" % \
                #    (n, basis, i, shell_nprim, atom_nprim)
                tst = ustart + atom_nprim
                tsp = ustart + atom_nprim + shell_nprim
                self.shells[shell_count] = GaussianShell(am, shell_nprim,
                    self.uoriginal_coefficients[tst:tsp],
                    self.ucoefficients[tst:tsp],
                    self.uerd_coefficients[tst:tsp],
                    self.uexponents[tst:tsp],
                    puream, n, xyz_ptr, bf_count)
                for thisbf in range(thisshell.nfunction()):
                    self.function_to_shell[bf_count] = shell_count
                    self.function_center[bf_count] = n
                    bf_count += 1
                for thisao in range(thisshell.ncartesian()):
                    self.ao_to_shell[ao_count] = shell_count
                    ao_count += 1
                atom_nprim += shell_nprim
                shell_count += 1

            # TODO huh?
            #xyz = self.molecule.xyz(n)
            #xyz_ptr[0] = xyz[0];
            #xyz_ptr[1] = xyz[1];
            #xyz_ptr[2] = xyz[2];
            #xyz_ptr += 3;

            if atom_nprim != uend - ustart:
                raise ValidationError("Problem with nprimitive in basis set construction!")

    @classmethod
    def build(cls, molecule, shells):
        """Builder factory method
        * @param molecule the molecule to build the BasisSet around
        * @param shells array of *atom-numbered* GaussianShells to build the BasisSet from
        * @return BasisSet corresponding to this molecule and set of shells

        """
        raise FeatureNotImplemented('BasisSet::build')
        #static boost::shared_ptr<BasisSet> build(boost::shared_ptr<Molecule> molecule, const std::vector<ShellInfo> &shells);
        #{
        #    //TODO fixme!!!
        #    boost::shared_ptr<BasisSet> basis(new BasisSet());
        #//    basis->molecule_ = molecule;
        #//    basis->shells_ = shells;
        #//    basis->refresh();
        #
        #    return basis;
        #}

    def constructor_basisset_center(self, bs, center):
        """
        * Creates a new basis set object for an atom, from an existing basis set
        * bs: the basis set to copy data from
        * center: the atom in bs to copy over

        """
        raise FeatureNotImplemented('BasisSet::constructor_basisset_center')
        #BasisSet::BasisSet(const BasisSet *bs, const int center)
        #{
        #    // Singletons; these should've been initialized by this point, but just in case
        #    if (initialized_shared_ == false)
        #        initialize_singletons();
        #    initialized_shared_ = true;
        #
        #    /*
        #     * First, find the shells we need, and grab the data
        #     */
        #    std::vector<double> uexps;
        #    std::vector<double> ucoefs;
        #    std::vector<double> uoriginal_coefs;
        #    std::vector<double> uerd_coefs;
        #    name_ = bs->name();
        #    n_shells_ = 0;
        #    n_uprimitive_ = 0;
        #    nao_ = 0;
        #    nbf_ = 0;
        #    for(int shelln = 0; shelln < bs->nshell(); ++shelln){
        #        const GaussianShell &shell = bs->shell(shelln);
        #        if(shell.ncenter() == center){
        #            int nprim = shell.nprimitive();
        #            for(int prim = 0; prim < nprim; ++prim){
        #                uexps.push_back(shell.exp(prim));
        #                ucoefs.push_back(shell.coef(prim));
        #                uoriginal_coefs.push_back(shell.original_coef(prim));
        #                uerd_coefs.push_back(shell.erd_coef(prim));
        #                n_uprimitive_++;
        #            }
        #            n_shells_++;
        #            nao_ += shell.ncartesian();
        #            nbf_ += shell.nfunction();
        #        }
        #    }
        #    nprimitive_ = n_uprimitive_;
        #
        #
        #    // Create a "molecule", i.e., an atom
        #    boost::shared_ptr<Molecule> mol = bs->molecule();
        #    molecule_ = boost::shared_ptr<Molecule>(new Molecule);
        #    int Z = mol->Z(center);
        #    double x = mol->x(center);
        #    double y = mol->y(center);
        #    double z = mol->z(center);
        #    double mass = mol->mass(center);
        #    double charge = mol->charge(center);
        #    std::string lab = mol->label(center);
        #    char* label = new char[lab.length() + 1];
        #    strcpy(label,lab.c_str());
        #    //Put the atomic info into mol
        #    molecule_->add_atom(Z, 0.0, 0.0, 0.0, label, mass, charge);
        #
        #
        #    /*
        #     * Allocate arrays
        #     */
        #    n_prim_per_shell_ = new int[n_shells_];
        #    // The unique primitives
        #    uexponents_ = new double[n_uprimitive_];
        #    ucoefficients_ = new double[n_uprimitive_];
        #    uoriginal_coefficients_ = new double[n_uprimitive_];
        #    uerd_coefficients_ = new double[n_uprimitive_];
        #    for(int i = 0; i < n_uprimitive_; ++i){
        #        uexponents_[i] = uexps[i];
        #        ucoefficients_[i] = ucoefs[i];
        #        uoriginal_coefficients_[i] = uoriginal_coefs[i];
        #        uerd_coefficients_[i] = uoriginal_coefs[i];
        #    }
        #
        #    shell_first_ao_ = new int[n_shells_];
        #    shell_first_basis_function_ = new int[n_shells_];
        #    shells_ = new GaussianShell[n_shells_];
        #    ao_to_shell_ = new int[nao_];
        #    function_to_shell_ = new int[nbf_];
        #    function_center_ = new int[nbf_];
        #    shell_center_ = new int[n_shells_];
        #    center_to_nshell_ = new int[1];
        #    center_to_shell_ = new int[1];
        #    xyz_ = new double[3];
        #
        #    /*
        #     * Now loop over shell for this atom, and point to the appropriate unique data
        #     */
        #    int shell_count = 0;
        #    int ao_count = 0;
        #    int bf_count = 0;
        #    puream_ = false;
        #    max_am_ = 0;
        #    max_nprimitive_ = 0;
        #    int prim_count = 0;
        #    for(int shelln = 0; shelln < bs->nshell(); ++shelln){
        #        const GaussianShell &shell = bs->shell(shelln);
        #        if(shell.ncenter() == center){
        #            center_to_nshell_[0] = n_shells_;
        #            center_to_shell_[0] = shell_count;
        #            shell_first_ao_[shell_count] = ao_count;
        #            shell_first_basis_function_[shell_count] = bf_count;
        #            int shell_nprim = shell.nprimitive();
        #            int am = shell.am();
        #            max_nprimitive_ = shell_nprim > max_nprimitive_ ? shell_nprim : max_nprimitive_;
        #            max_am_ = max_am_ > am ? max_am_ : am;
        #            shell_center_[shell_count] = center;
        #            GaussianType puream = shell.is_pure() ? Pure : Cartesian;
        #            if(puream)
        #                puream_ = true;
        #            shells_[shell_count] = GaussianShell(am, shell_nprim, &uoriginal_coefficients_[prim_count],
        #                    &ucoefficients_[prim_count], &uerd_coefficients_[prim_count], &uexponents_[prim_count], puream, center, xyz_, bf_count);
        #            for(int thisbf = 0; thisbf < shell.nfunction(); ++thisbf){
        #                function_to_shell_[bf_count] = shell_count;
        #                function_center_[bf_count++] = center;
        #            }
        #            for(int thisao = 0; thisao < shell.ncartesian(); ++thisao){
        #                ao_to_shell_[ao_count++] = shell_count;
        #            }
        #            shell_count++;
        #            prim_count += shell_nprim;
        #        }
        #    }
        #    xyz_[0] = xyz_[1] = xyz_[2] = 0.0;
        #}

    def initialize_singletons(self):
        """Initialize singleton values that are shared by all basis set objects."""
        # Populate the exp_ao arrays
        for l in range(self.LIBINT_MAX_AM):
            for i in range(l + 1):
                x = l - i
                for j in range(i + 1):
                    y = i - j
                    z = j
                    self.exp_ao[l].append([x, y, z])

    # <<< Simple Methods for Coordinates >>>
    # <<< Simple Methods for Fragmentation >>>
    # <<< Methods for Construction >>>

    # <<< Simple Methods for Basic BasisSet Information >>>

    def name(self):
        """Returns the name of this basis set"""
        return self.name

    def set_name(self, name):
        """Sets the name of this basis set"""
        self.name = name

    def nprimitive(self):
        """Number of primitives.
        *  @return The total number of primitives in all contractions.

        """
        return self.PYnprimitive

    def max_nprimitive(self):
        """Maximum number of primitives in a shell.
        *  Examines each shell and find the shell with the maximum number of primitives returns that
        *  number of primitives.
        *  @return Maximum number of primitives.

        """
        return self.max_nprimitive

    def nshell(self):
        """Number of shells.
        *  @return Number of shells.

        """
        return self.n_shells

    def nao(self):
        """Number of atomic orbitals (Cartesian).
        * @return The number of atomic orbitals (Cartesian orbitals, always).

        """
        return self.PYnao

    def nbf(self):
        """Number of basis functions (Spherical).
        *  @return The number of basis functions (Spherical, if has_puream() == true).

        """
        return self.PYnbf

    def max_am(self):
        """Maximum angular momentum used in the basis set.
        *  @return Maximum angular momentum.

        """
        return self.PYmax_am

    def has_puream(self):
        """Spherical harmonics?
        *  @return true if using spherical harmonics

        """
        return self.puream

    def max_function_per_shell(self):
        """Compute the maximum number of basis functions contained in a shell.
        *  @return The max number of basis functions in a shell.

        """
        return 2 * self.PYmax_am + 1 if self.puream else (self.PYmax_am + 1) * (self.PYmax_am + 2) / 2

    def molecule(self):
        """Molecule this basis is for.
        *  @return Shared pointer to the molecule for this basis set.

        """
        return self.molecule

    def shell_to_ao_function(self, i):
        """Given a shell what is its first AO function
        *  @param i Shell number
        *  @return The function number for the first function for the i'th shell.

        """
        return self.shell_first_ao[i]

    def shell_to_center(self, i):
        """Given a shell what is its atomic center
        *  @param i Shell number
        *  @return The atomic center for the i'th shell.

        """
        return self.shell_center[i]

    def shell_to_basis_function(self, i):
        """Given a shell what is its first basis function (spherical) function
        *  @param i Shell number
        *  @return The function number for the first function for the i'th shell.

        """
        return self.shell_first_basis_function[i]

    def function_to_shell(self, i):
        """Given a function number what shell does it correspond to."""
        return self.function_to_shell[i]

    def function_to_center(self, i):
        """Given a function what is its atomic center
        *  @param i Function number
        *  @return The atomic center for the i'th function.

        """
        return self.function_center[i]

    def ao_to_shell(self, i):
        """Given a Cartesian function (AO) number what shell does it correspond to."""
        return self.ao_to_shell[i]

    def shell(self, si, center=None):
        """Return the si'th Gaussian shell on center
        *  @param i Shell number
        *  @return A shared pointer to the GaussianShell object for the i'th shell.

        """
        if center is not None:
            si += self.center_to_shell[center]
        if si < 0 or si > self.nshell():
            text = """BasisSet::shell(si = %d), requested a shell out-of-bound.\n   Max shell size: %d\n   Name: %s\n""" % \
                (si, self.nshell(), self.name())
            raise ValidationError("BasisSet::shell: requested shell is out-of-bounds:\n%s" % (text))
        return self.shells[si]

    # <<< Methods for Printing >>>

    def print_by_level(self, out=None, level=2):
        """Print basis set information according to the level of detail in print_level
        *  @param out The file stream to use for printing. Defaults to outfile.
        *  @param print_level: < 1: Nothing
                                 1: Brief summary
                                 2: Summary and contraction details
                               > 2: Full details
                               Defaults to 2

        """
        if level < 1:
            return
        elif level == 1:
            text = self.pyprint(out=None)
        elif level == 2:
            text = self.print_summary(out=None)
        elif level > 2:
            text = self.print_detail(out=None)

        if out is None:
            print text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    def pyprint(self, out=None):
        """Print the basis set.
        *  @param out The file stream to use for printing. Defaults to outfile.

        """
        text = ''
        text += """  Basis Set: %s\n""" % (self.name)
        text += """    Number of shells: %d\n""" % (self.nshell())
        text += """    Number of basis function: %d\n""" % (self.nbf())
        text += """    Number of Cartesian functions: %d\n""" % (self.nao())
        text += """    Spherical Harmonics?: %s\n""" % ('true' if self.has_puream() else 'false')
        text += """    Max angular momentum: %d\n\n""" % (self.max_am())

        if out is None:
            return text
        else:
            with open(outfile, mode='w') as handle:
                handle.write(text)

    def print_summary(self, out=None):
        """Prints a short string summarizing the basis set
        *  @param out The file stream to use for printing. Defaults to outfile.

        """
        text = ''
        text += """  -AO BASIS SET INFORMATION:\n"""
        text += """    Name                   = %s\n""" % (self.name)
        text += """    Total number of shells = %d\n""" % (self.nshell())
        text += """    Number of primitives   = %d\n""" % (self.nprimitive())
        text += """    Number of AO           = %d\n""" % (self.nao())
        text += """    Number of SO           = %d\n""" % (self.nbf())
        text += """    Maximum AM             = %d\n""" % (self.max_am())
        text += """    Spherical Harmonics    = %s\n""" % ('TRUE' if self.puream else 'FALSE')
        text += """\n"""
        text += """  -Contraction Scheme:\n"""
        text += """    Atom   Type   All Primitives // Shells:\n"""
        text += """   ------ ------ --------------------------\n"""

        print self.molecule.natom()
        self.molecule.print_out()
        for A in range(self.molecule.natom()):

            nprims = [0] * (self.PYmax_am + 1)
            nunique = [0] * (self.PYmax_am + 1)
            nshells = [0] * (self.PYmax_am + 1)
            amtypes = [None] * (self.PYmax_am + 1)

            text += """    %4d    """ % (A + 1)
            text += """%2s     """ % (self.molecule.symbol(A))

            first_shell = self.center_to_shell[A]
            n_shell = self.center_to_nshell[A]

            for Q in range(n_shell):
                shell = self.shells[Q + first_shell]
                nshells[shell.am()] += 1
                nunique[shell.am()] += shell.nprimitive()
                nprims[shell.am()] += shell.nprimitive()
                amtypes[shell.am()] = shell.amchar()

            # All Primitives
            for l in range(self.PYmax_am + 1):
                if nprims[l] == 0:
                    continue
                text += """%d%c """ % (nprims[l], amtypes[l])

            # Shells
            text += """// """
            for l in range(self.PYmax_am + 1):
                if nshells[l] == 0:
                    continue
                text += """%d%c """ % (nshells[l], amtypes[l])
            text += """\n"""
        text += """\n"""

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    def print_detail(self, out=None):
        """Prints a detailed PSI3-style summary of the basis (per-atom)
        *  @param out The file stream to use for printing. Defaults to outfile.

        """
        text = self.print_summary(out=None)

        text += """  ==> AO Basis Functions <==\n"""
        text += '\n'
        text += """    [ %s ]\n""" % (self.name)
        text += """    spherical\n""" if self.has_puream() else """    cartesian\n"""
        text += """    ****\n"""

        for uA in range(self.molecule.nunique()):
            A = self.molecule.unique(uA)
            text += """   %2s %3d\n""" % (self.molecule.symbol(A), A + 1)
            first_shell = self.center_to_shell[A]
            n_shell = self.center_to_nshell[A]

            for Q in range(n_shell):
                text += self.shells[Q + first_shell].pyprint(outfile=None)
            text += """    ****\n"""
        text += """\n"""

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    def print_detail_cfour(self, out=None):
        """Returns a string in CFOUR-style of the basis (per-atom)
        *  Format from http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.OldFormatOfAnEntryInTheGENBASFile

        """
        text = ''

        for uA in range(self.molecule.nunique()):
            A = self.molecule.unique(uA)
            text += """%s:P4_%d\n""" % (self.molecule.symbol(A), A + 1)
            text += """PSI4 basis %s for element %s atom %d\n\n""" % \
                (self.name.upper(), self.molecule.symbol(A), A + 1)

            first_shell = self.center_to_shell[A]
            n_shell = self.center_to_nshell[A]

            max_am_center = 0
            for Q in range(n_shell):
                max_am_center = self.shells[Q + first_shell].am() if \
                self.shells[Q + first_shell].am() > max_am_center else max_am_center

            shell_per_am = [[] for i in range(max_am_center + 1)]
            for Q in range(n_shell):
                shell_per_am[self.shells[Q + first_shell].am()].append(Q)

            # Write number of shells in the basis set
            text += """%3d\n""" % (max_am_center + 1)

            # Write angular momentum for each shell
            for am in range(max_am_center + 1):
                text += """%5d""" % (am)
            text += '\n'

            # Write number of contracted basis functions for each shell
            for am in range(max_am_center + 1):
                text += """%5d""" % (len(shell_per_am[am]))
            text += '\n'

            exp_per_am = [[] for i in range(max_am_center + 1)]
            coef_per_am = [[] for i in range(max_am_center + 1)]
            for am in range(max_am_center + 1):
                # Collect unique exponents among all functions
                for Q in range(len(shell_per_am[am])):
                    for K in range(self.shells[shell_per_am[am][Q] + first_shell].nprimitive()):
                        if self.shells[shell_per_am[am][Q] + first_shell].exp(K) not in exp_per_am[am]:
                            exp_per_am[am].append(self.shells[shell_per_am[am][Q] + first_shell].exp(K))

                # Collect coefficients for each exp among all functions, zero otherwise
                for Q in range(len(shell_per_am[am])):
                    K = 0
                    for ep in range(len(exp_per_am[am])):
                        if abs(exp_per_am[am][ep] - self.shells[shell_per_am[am][Q] + first_shell].exp(K)) < 1.0e-8:
                            coef_per_am[am].append(self.shells[shell_per_am[am][Q] + first_shell].original_coef(K))
                            if (K + 1) != self.shells[shell_per_am[am][Q] + first_shell].nprimitive():
                                K += 1
                        else:
                            coef_per_am[am].append(0.0)

            # Write number of exponents for each shell
            for am in range(max_am_center + 1):
                text += """%5d""" % (len(exp_per_am[am]))
            text += '\n\n'

            for am in range(max_am_center + 1):
                # Write exponents for each shell
                for ep in range(len(exp_per_am[am])):
                    text += """%14.7f""" % (exp_per_am[am][ep])
                    if ((ep + 1) % 5 == 0) or ((ep + 1) == len(exp_per_am[am])):
                        text += '\n'
                text += '\n'

                # Write contraction coefficients for each shell
                for ep in range(len(exp_per_am[am])):
                    for bf in range(len(shell_per_am[am])):
                        text += """%10.7f """ % (coef_per_am[am][bf * len(exp_per_am[am]) + ep])
                    text += '\n'
                text += '\n'

        if out is None:
            return text
        else:
            with open(out, mode='w') as handle:
                handle.write(text)

    def refresh(self):
        """Refresh internal basis set data. Useful if someone has pushed
        to shells_. Pushing to shells_ happens in the BasisSetParsers, so
        the parsers will call refresh(). This function is now defunct.

        """
        raise FeatureNotImplemented('BasisSet::refresh')

    def nshell_on_center(self, i):
        """Return the number of shells on a given center."""
        return self.center_to_nshell[i]

    def shell_on_center(self, center, shell):
        """Return the overall shell number"""
        return self.center_to_shell[center] + shell

    def atomic_basis_set(self, center):
        """Return a BasisSet object containing all shells at center i
        * Used for Atomic HF computations for SAD Guesses
        * @param center Atomic center to provide a basis object for.
        * @returns A new basis set object for the atomic center.

        """
        raise FeatureNotImplemented('BasisSet::atomic_basis_set')

    @staticmethod
    def zero_ao_basis_set(cls):
        """Returns an empty basis set object.
        Returns a BasisSet object that actually has a single s-function
        at the origin with an exponent of 0.0 and contraction of 1.0.
        *  @return A new empty BasisSet object.

        """
        # TODO check
        # In the new implementation, we simply call the default constructor
        return BasisSet()

    def zero_so_basis_set(cls, factory):
        """ **NYI** Returns an empty SO basis set object.
        *  Returns an SOBasis object that actually has a single s-function
        *  at the origin with an exponent of 0.0 and contraction of 1.0.
        *  @return A new empty SOBasis object.

        """
        raise FeatureNotImplemented('BasisSet::zero_so_basis_set')  # FINAL

    @classmethod
    def test_basis_set(cls, max_am):
        """Returns a shell-labeled test basis set object
        * @param max_am maximum angular momentum to build
        * @return pair containing shell labels and four-center
        * test basis for use in benchmarking
        * See libmints/benchmark.cc for details

        """
        raise FeatureNotImplemented('BasisSet::test_basis_set')
        #max_centers = 4
        #max_primitives = 10

        #nprim = [10, 1, 6, 1, 2, 1, 1, 1, 1, 1]
        #am = [0, 0, 1, 1, 2, 2, 3, 4, 5, 2]

        #c = []
        #c.append([0.0] * 10)
        #c.append([0.0])
        #c.append([0.0] * 6)
        #c.append([0.0])
        #c.append([0.0] * 2)
        #c.append([0.0])
        #c.append([0.0])
        #c.append([0.0])
        #c.append([0.0])
        #c.append([0.0])

        #e = []
        #e.append([0.0] * 10)
        #e.append([0.0])
        #e.append([0.0] * 6)
        #e.append([0.0])
        #e.append([0.0] * 2)
        #e.append([0.0])
        #e.append([0.0])
        #e.append([0.0])
        #e.append([0.0])
        #e.append([0.0])

        #c[0][0] = 0.458878E-03
        #c[0][1] = 0.355070E-02
        #c[0][2] = 0.182618E-01
        #c[0][3] = 0.716650E-01
        #c[0][4] = 0.212346E+00
        #c[0][5] = 0.416203E+00
        #c[0][6] = 0.373020E+00
        #c[0][7] = 0.625054E-01
        #c[0][8] = 0.624532E-02
        #c[0][9] = 0.243374E-02
        #c[1][0] =     1.0
        #c[2][0] = 0.458878E-03
        #c[2][1] = 0.355070E-02
        #c[2][2] = 0.182618E-01
        #c[2][3] = 0.716650E-01
        #c[2][4] = 0.212346E+00
        #c[2][5] = 0.416203E+00
        #c[3][0] =     1.0
        #c[4][0] = 0.458878E-03
        #c[4][1] = 0.355070E-02
        #c[5][0] =     1.0
        #c[6][0] =     1.0
        #c[7][0] =     1.0
        #c[8][0] =     1.0
        #c[9][0] =     1.0

        #e[0][0] = 31700.0
        #e[0][1] =  4755.0
        #e[0][2] =  1082.0
        #e[0][3] =   306.0
        #e[0][4] =    99.0
        #e[0][5] =    33.0
        #e[0][6] =    13.0
        #e[0][7] =     4.0
        #e[0][8] =     2.0
        #e[0][9] =     0.5
        #e[1][0] =     1.0
        #e[2][0] = 31700.0
        #e[2][1] =  4755.0
        #e[2][2] =  1082.0
        #e[2][3] =   306.0
        #e[2][4] =    99.0
        #e[2][5] =    33.0
        #e[3][0] =     1.0
        #e[4][0] = 31700.0
        #e[4][1] =  4755.0
        #e[5][0] =     1.0
        #e[6][0] =     1.0
        #e[7][0] =     1.0
        #e[8][0] =     1.0
        #e[9][0] =     1.0

        #labels = []
        #if max_am > -1:
        #    labels.append("S")
        #    labels.append("s")
        #    max_shells = 2
        #if max_am > 0:
        #    labels.append("P")
        #    labels.append("p")
        #    max_shells = 4
        #if max_am > 1:
        #    labels.append("D")
        #    labels.append("d")
        #    max_shells = 6
        #if max_am > 2:
        #    labels.append("f")
        #    max_shells = 7
        #if max_am > 3:
        #    labels.append("g")
        #    max_shells = 8
        #if max_am > 4:
        #    labels.append("h")
        #    max_shells = 9
        #if max_am > 5:
        #    labels.append("i")
        #    max_shells = 10

        #new_basis = BasisSet()

        ## Add 4 atoms to the molecule for this basis (max integal centers is 4 at the moment)
        #new_basis.molecule = Molecule()
        ## Ghost atoms are now handled differently, they are not added to the normal xyz information array,
        ## but to the fxyz array.
        #x = 0.0
        #for A in range(max_centers):
        #    new_basis.molecule.add_atom(0, x, x, x)
        #    x += 1.0

        ## Setup all the parameters needed for a zero basis set
        #new_basis.shell_center = [0] * max_shells * max_centers
        #for A in range(max_centers):
        #    for Q in range(max_shells):
        #        new_basis.shell_center[A * max_shells + Q] = A
        #new_basis.max_nprimitive = max_primitives
        #new_basis.PYmax_am = max_am

        ## We'll time puream for now
        #new_basis.puream = True

        ## Add shells
        #for A in range(max_centers):
        #    center = new_basis.molecule.fxyz(A)
        #    for Q in range(max_shells):
        #        new_basis.shells.append(GaussianShell(am[Q], nprim[Q], c[Q], e[Q], 'Pure', A, center, 0))
        #        #ShellInfo(self, am, c, e, pure, nc, center, start, pt='Normalized'):
        #new_basis.refresh()

        #print labels
        #print new_basis
        ##return make_pair(labels, new_basis)
        #return new_basis

#        self.name = None
##        self.shells = None
##        self.molecule = None
#        self.PYnao = None
#        self.PYnbf = None
#        self.n_uprimitive = None
#        self.n_shells = None
#        self.PYnprimitive = None
##        self.PYmax_am = None
##        self.max_nprimitive = None
##        self.puream = None
#        self.n_prim_per_shell = None
#        self.shell_first_ao = None
#        self.shell_first_basis_function = None
##        self.shell_center = None
#        self.function_to_shell = None
#        self.ao_to_shell = None
#        self.function_center = None
#        self.center_to_nshell = None
#        self.center_to_shell = None
#        self.uexponents = None
#        self.ucoefficients = None
#        self.uoriginal_coefficients = None
#        self.uerd_coefficients = None
#        self.xyz = None

    @classmethod
    def construct(cls, parser, mol, role):
        """Returns a new BasisSet object.
        * Returns a new BasisSet object configured with the provided Molecule object.
        * @param parser The basis set parser object that will be used to interpret the basis set file.
        * @param mol Molecule to construct basis set for.
        * @param basisnames Name of the basis set for each atom in molecule to search for in pbasis.dat
        * @return A new basis set object constructed from the information passed in.

        """
        # Update geometry in molecule, if there is a problem an exception is thrown.
        mol.update_geometry()

        # Paths to search for gbs files
        basisPath = ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':')]) + \
            ':' + os.path.abspath(os.environ.get('PSIDATADIR')) + '/basis'
            # TODO: submission dir, #psi4.Process.environment["PSIDATADIR"] + '/basis' + ':' + psi4.psi_top_srcdir() + '/lib/basis'
            #std::string psiPath = PSIOManager::shared_object()->get_default_path() +
            # somehow the calling dir is included, i'm not seeing how

        # Map of GaussianShells
        basis_atom_shell = collections.OrderedDict()  # basis_atom_shell[str(basis)][str(atom)] = GaussianShells(shells)
        names = collections.OrderedDict()  # names[str(basis)] = int(found)

        for atom in range(mol.natom()):
            symbol = mol.atom_entry(atom).symbol()  # O, He
            label = mol.atom_entry(atom).label()  # O3, C_Drot, He
            try:
                basisname = mol.atom_entry(atom).basisset(role)
            except KeyError:
                raise BasisSetNotDefined("""BasisSet::construct: No basis set specified for %s and %s.""" %
                    (symbol, role))

            names[basisname] = 1
            # Add basisname, symbol to the list by clearing the vector.
            if basisname not in basis_atom_shell:
                basis_atom_shell[basisname] = collections.OrderedDict()
            basis_atom_shell[basisname][label] = []

        for basisfirst, basissecond in basis_atom_shell.items():
            filename = cls.make_filename(basisfirst)
            fullfilename = search_file(filename, basisPath)

            if fullfilename is None:
                print '\nBasis set file %s failed to load\n\n' % (filename)
                print '\nSearch path that was tried:\n'
                print ', '.join(map(str, basisPath.split(':')))
                raise BasisSetFileNotFound('Basis set file loading problem for ' + (filename))

            lines = parser.load_file(fullfilename)
            for atomfirst, atomsecond in basissecond.items():
                label = atomfirst
                symbol = re.split('\d|_', label)[0]

                # Tries to find basis for label (e.g., N88) in basis. Failing that,
                #   tries to find basis for symbol (e.g., N) in basis.
                try:
                    basis_atom_shell[basisfirst][label] = parser.parse(label, lines)
                except BasisSetNotFound:
                    if label == symbol:
                        raise BasisSetNotFound
                    else:
                        basis_atom_shell[basisfirst][label] = parser.parse(symbol, lines)

        basisset = BasisSet(role, mol, basis_atom_shell)

        basisset.name = ''
        for name in names:
            basisset.name += name + ' + '
        if basisset.name.endswith(' + '):
            basisset.name = basisset.name[:-3]

        return basisset

    @staticmethod
    def make_filename(name):
        """Converts basis set name to a compatible filename.
        * @param basisname Basis name
        * @return Compatible file name.

        """
        # Modify the name of the basis set to generate a filename: STO-3G -> sto-3g
        basisname = name

        # First make it lower case
        basisname = basisname.lower()

        # Replace all '(' with '_'
        basisname = basisname.replace('(', '_')

        # Replace all ')' with '_'
        basisname = basisname.replace(')', '_')

        # Replace all ',' with '_'
        basisname = basisname.replace(',', '_')

        # Replace all '*' with 's'
        basisname = basisname.replace('*', 's')

        # Replace all '+' with 'p'
        basisname = basisname.replace('+', 'p')

        # Add file extension
        basisname += '.gbs'

        return basisname

    def get_ao_sorted_shell(self, i):
        """Returns the value of the sorted shell list. Defunct"""
        raise FeatureNotImplemented('BasisSet::get_ao_sorted_shell')

    def get_ao_sorted_list(self):
        """Returns the vector of sorted shell list. Defunct"""
        raise FeatureNotImplemented('BasisSet::get_ao_sorted_list')

    def compute_phi(self, phi_ao, x, y, z):
        """Returns the values of the basis functions at a point"""

        phi_ao = [0.0] * self.nao()
        ao = 0
        for ns in range(self.nshell()):
            shell = self.shells[ns]
            am = shell.am()
            nprim = shell.nprimitive()
            a = shell.exps()
            c = shell.coefs()

            xyz = shell.center()
            dx = x - xyz[0]
            dy = y - xyz[1]
            dz = z - xyz[2]
            rr = dx * dx + dy * dy + dz * dz

            cexpr = 0
            for np in range(nprim):
                cexpr += c[np] * math.exp(-a[np] * rr)

            for l in range(INT_NCART(am)):
                components = exp_ao[am][l]
                phi_ao[ao + l] += pow(dx, components[0]) * \
                                pow(dy, components[1]) * \
                                pow(dz, components[2]) * \
                                cexpr

            ao += INT_NCART(am)

#    /** Concatenates two basis sets together into a new basis without reordering anything.
#     *  Unless you know what you're doing, you should use the '+' operator instead of
#     *  this method.
#     */
#    BasisSet concatenate(const BasisSet& b) const;

#    boost::shared_ptr<BasisSet> concatenate(const boost::shared_ptr<BasisSet>& b) const;

#    /** Concatenates two basis sets together into a new basis without reordering anything.
#     *  Unless you know what you're doing, you should use the '+' operator instead of
#     *  this method.
#     */
#    //static boost::shared_ptr<BasisSet> concatenate(const boost::shared_ptr<BasisSet>& a, const boost::shared_ptr<BasisSet>& b);

#    /** Adds this plus another basis set and returns the result. Equivalent to the '+' operator.
#     */
#    BasisSet add(const BasisSet& b) const;

#    boost::shared_ptr<BasisSet> add(const boost::shared_ptr<BasisSet>& b) const;

#    // BasisSet friends
#    friend class Gaussian94BasisSetParser;
#    friend BasisSet operator +(const BasisSet& a, const BasisSet& b);
#    friend boost::shared_ptr<BasisSet> operator +(const boost::shared_ptr<BasisSet>& a, const boost::shared_ptr<BasisSet>& b);

#    // Adds 2 shared basis set objects together
#    static boost::shared_ptr<BasisSet> add(const boost::shared_ptr<BasisSet>& a, const boost::shared_ptr<BasisSet>& b) {
#        return boost::shared_ptr<BasisSet>(new BasisSet(*a.get() + *b.get()));
#    }


def shell_sorter_ncenter(d1, d2):
    #GaussianShell& d1, const GaussianShell& d2)
    return d1.ncenter() < d2.ncenter()


def shell_sorter_am(d1, d2):
    #const GaussianShell& d1, const GaussianShell& d2)
    return d1.am() < d2.am()

#BasisSet operator +(const BasisSet& a, const BasisSet& b) {
#    if (a.molecule() != b.molecule()) {
#        fprintf(stderr, "BasisSet::operator+ : Unable to add basis sets from different molecules.");
#        return BasisSet();
#    }
#    BasisSet temp;
#
#if 0 //TODO fixme!
#    temp.name_ = a.name_ + " + " + b.name_;
#    temp.molecule_ = a.molecule();
#
#    // Copy a's shells to temp
#    temp.shells_ = a.shells_;
#
#    // Append b's shells to temp
#    temp.shells_.insert(temp.shells_.end(), b.shells_.begin(), b.shells_.end());
#
#    // Sort by center number
#    std::sort(temp.shells_.begin(), temp.shells_.end(), shell_sorter_ncenter);
#
#    // Call refresh to regenerate center_to_shell and center_to_nshell
#    temp.refresh();
#
#    // Sort by AM in each center
#//    for (int atom=0; atom < temp.molecule_->natom(); ++atom) {
#//        std::sort(temp.shells_.begin()+temp.center_to_shell_[atom],
#//                  temp.shells_.begin()+temp.center_to_shell_[atom]+temp.center_to_nshell_[atom],
#//                  shell_sorter_am);
#//    }
#
#//    temp.refresh();
##endif
#    return temp;
#}

#inline
#boost::shared_ptr<BasisSet> operator +(const boost::shared_ptr<BasisSet>& a, const boost::shared_ptr<BasisSet>& b) {
#    return boost::shared_ptr<BasisSet>(new BasisSet(*a.get() + *b.get()));
#}
#}


#
#
#//boost::shared_ptr<SOBasisSet> BasisSet::zero_so_basis_set(const boost::shared_ptr<IntegralFactory>& factory)
#//{
#//    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
#//    boost::shared_ptr<SOBasisSet> sozero(new SOBasisSet(zero, factory));
#//    return sozero;
#//}
#
#
#
#
#boost::shared_ptr<BasisSet> BasisSet::atomic_basis_set(int center)
#{
#    return boost::shared_ptr<BasisSet>(new BasisSet(this, center));
#}
#
#
#BasisSet BasisSet::concatenate(const BasisSet& b) const {
#
#    BasisSet temp;
#//TODO fixme!!!
##if 0
#    temp.name_ = name_ + " + " + b.name_;
#    temp.molecule_ = molecule();
#
#    // Copy a's shells to temp
#    temp.shells_ = shells_;
#
#    // Append b's shells to temp
#    temp.shells_.insert(temp.shells_.end(), b.shells_.begin(), b.shells_.end());
#
#    // Call refresh to regenerate center_to_shell and center_to_nshell
#    temp.refresh();
##endif
#    return temp;
#}
#
#boost::shared_ptr<BasisSet> BasisSet::concatenate(const boost::shared_ptr<BasisSet>& b) const {
#    return boost::shared_ptr<BasisSet>(new BasisSet(concatenate(*b.get())));
#}
#
#BasisSet BasisSet::add(const BasisSet& b) const {
#
#    BasisSet temp;
#//TODO fixme!!
##if 0
#    temp.name_ = name_ + " + " + b.name_;
#    temp.molecule_ = molecule();
#
#    // Copy a's shells to temp
#//    temp.shells_ = shells_;
#
#    // Append b's shells to temp
#//    std::vector<GaussianShell>::const_iterator iter = b.shells_.begin();
#//    for (; iter != b.shells_.end(); iter++)
#//        temp.shells_.push_back(*iter);
#
#//    temp.shells_.insert(temp.shells_.end(), b.shells_.begin(), b.shells_.end());
#
#    // Loop over atoms
#    for (int atom=0; atom<molecule()->natom(); ++atom) {
#        for (int shella=0; shella<nshell_on_center(atom); ++shella)
#            temp.shells_.push_back(shell(atom, shella));
#
#        for (int shellb=0; shellb<b.nshell_on_center(atom); ++shellb)
#            temp.shells_.push_back(b.shell(atom, shellb));
#    }
#
#    // Call refresh to regenerate center_to_shell and center_to_nshell
#    temp.refresh();
#
#    // Sort by center number
#//    std::sort(temp.shells_.begin(), temp.shells_.end(), shell_sorter_ncenter);
#
#    // Call refresh to regenerate center_to_shell and center_to_nshell
#//    temp.refresh();
#
#    // Sort by AM in each center
#//    for (int atom=0; atom < temp.molecule_->natom(); ++atom) {
#//        std::sort(temp.shells_.begin()+temp.center_to_shell_[atom],
#//                  temp.shells_.begin()+temp.center_to_shell_[atom]+temp.center_to_nshell_[atom],
#//                  shell_sorter_am);
#//    }
##endif
#    return temp;
#}
#
#boost::shared_ptr<BasisSet> BasisSet::add(const boost::shared_ptr<BasisSet>& b) const {
#    return boost::shared_ptr<BasisSet>(new BasisSet(add(*b.get())));
#}
