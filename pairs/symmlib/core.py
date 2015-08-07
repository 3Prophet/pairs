import sympy
from sympy import abc
import re
from cctbx import sgtbx
import sys
import collections
import numpy
import math

class SpaceGroup(sgtbx.space_group):
    """Facade class creating simpler interface for sgtbx.space_group class"""
    def __init__(self, sg_nr):
        """Initialize with SG number"""
        self._nr = sg_nr
        self._sg_hall = sgtbx.space_group_symbols(sg_nr).hall()
        super(SpaceGroup, self).__init__(self._sg_hall)
        self._order = self.order_z()
        self._symops = []

    def __iter__(self):
        return iter(self.symops)
    
    def __repr__(self):
        return "SpaceGroup(sg_nr = {})".format(self.nr)

    @property
    def nr(self):
        return self._nr

    @property
    def order(self):
        return self._order

    @property
    def symops(self):
        """Symmetry operations in symbolic form"""
        if not self._symops:
            for i in range(self.order):
                self._symops.append(self(i).as_xyz())
        return self._symops

class Parser(list):
    """Example of what it does: '-x,2y,z+1/2'--->[['-1','x','0'],
                                                  ['2','y','0'],
                                                  ['1','z','1/2']]
    It returns list which entries contain results of parsing.
    """

    def __init__(self, symop):
        """Symop has to be of the form '-x,2y,z' """
        self._symop = symop
        self._parsed = self.__call__()
        super(Parser, self).__init__(self._parsed)
        
    def __call__(self):
        """Returns the final result upon the call"""
        tags = self._tags()
        rt_pairs_start = self._rt_pairs_start(tags)
        rt_mult_order = self._rt_mult_order(rt_pairs_start)
        return [item for item in rt_mult_order]

    def _tags(self):
        """  1-x,y+1,-z'--> ['1-x','y+1','-z']"""
        for item in self._symop.strip().split(','):
            yield item

    def _rt_pairs_start(self, tags):
        """ '-1/2x+1' --->['-1/2x','+1']
                 '-x' --->['-x']
                 etc...
        """
        for tag in tags:
            rt_list_start = re.split(r'([+-]?(?:\d+)?/?(?:\d+)?[xyz])',
                                        tag.strip())
            cleaned = filter(lambda item: item != "", rt_list_start)
            yield cleaned

    def _rt_mult_order(self, rt_list_start):
        """ ['-1/2x','+1'] ---> ['-1/2','x','+1'] 
                    ['-x'] ---> ['-1','x','0']
                                [multiplier, [xyz], translation]   
        """
        test_set = set('xyz')
        for alist in rt_list_start:
            rt_mult_order = [None, None, None] #[multiplier, [xyz], translation]
            for item in alist:
                if not test_set.isdisjoint(set(item)):  #if contains x,y or z 
                    patt = re.compile(r'(?P<sign>[+-]?)(?P<mult>[+-]?(?:\d+)?/?(?:\d+)?)(?P<let>[xyz])')
                    parsed = re.search(patt, item)
                    if parsed.group('mult') == '':
                        if parsed.group('sign') == '':
                            rt_mult_order[0] = '1'
                        else:
                            rt_mult_order[0] = parsed.group('sign')+'1'
                    else:
                        if parsed.group('sign') == '':
                            rt_mult_order[0] = parsed.group('mult')
                        else:
                            rt_mult_order[0] = parsed.group('sign')+\
                                            parsed.group('mult')
                    rt_mult_order[1] = parsed.group('let')
                else:                                   #if translation
                    rt_mult_order[2] = item
            if rt_mult_order[-1] == None:               #if no translation
                rt_mult_order[-1]  = '0'
            yield rt_mult_order

    @classmethod
    def tests(Cls):
        """Some tests, fill free to fill sym.operations of interest into
        symops for testing."""
        symops = ['x,y,z', 'x,z,y','-x,y,z', '-2x,y,z', 'x,y,1/2-z',
                'x,y,1/2-2z', '1/2x,y,z']
        print("\nRUNNING TESTS FOR CLASS {}...\n".format(Cls.__name__))
        for symop in symops:
            print("Parsing {} |--->".format(symop)), Cls(symop)
        print("\nDONE!\n")

class MetricTensor(object):
    """Metric tensor object for calculating lengths, angles and distances"""
    def __init__(self, lattice = (1.,1.,1.,90.,90.,90.), tensor = None):
        if tensor != None:
            if not isinstance(tensor, numpy.ndarray):
                raise TypeError("Tensor should be an instance of numpy.ndarray !")
            if tensor.shape != (3,3):
                raise IndexError("Shape of the tensor should be (3,3)")
            self._tensor = tensor
        else:

            if not isinstance(lattice, tuple):
                raise TypeError("Lattice has to be a tuple !")
            if not len(lattice) == 6:
                raise IndexError("Lattice should contain 6 entries !")
            a,b,c,alpha,beta,gamma = lattice
            cos_alpha = math.cos(math.radians(alpha))
            cos_beta = math.cos(math.radians(beta))
            cos_gamma = math.cos(math.radians(gamma))
            self._lattice = lattice
            self._tensor = numpy.array([
                                [a*a,           a*b*cos_gamma, a*c*cos_beta],
                                [b*a*cos_gamma, b*b,           b*c*cos_alpha],
                                [c*a*cos_beta,  c*b*cos_alpha, c*c]
                                    ])

    @property 
    def tensor(self):
        return self._tensor
    @tensor.setter
    def tensor(self, value):
        raise NotImplementedError
    @tensor.deleter
    def tensor(self, value):
        raise NotImplementedError

    def length(self, vec):
        """Calculates length of the vector"""
        vec_num = eval(self.__class__.__name__)()._in2numpyvec(vec)
        length =numpy.dot(vec_num,numpy.dot(self.tensor, vec_num.T))[0][0]**0.5
        return length

    def distance(self, vec1, vec2):
        """Returns distance between 2 points given by vec1 and vec2"""
        vec_num1 = eval(self.__class__.__name__)()._in2numpyvec(vec1)
        vec_num2 = eval(self.__class__.__name__)()._in2numpyvec(vec2)
        diff = vec_num2 - vec_num1
        distance = numpy.dot(diff, numpy.dot(self.tensor, diff.T))[0][0]**0.5
        return distance

    def angle(self, vec1, vec2):
        raise NotImplementedError

    def inverted(self):
        """Returns an instance of the inverted tensor"""
        return self.__class__(tensor = numpy.linalg.inv(self.tensor))


    @staticmethod
    def _in2numpyvec(vec):
        """Transforms vector in form of tuple, or set, or list into numpy 
        vector of shape: (1,3).
        Examples of vec: [1,2,3] or (1,2,3)- i.e. vec should be nonnested 
        container of lengtn 3."""
        if  not any([isinstance(vec, list),  isinstance(vec, tuple),
            isinstance(vec, set)]):
            raise TypeError("Allowed types are: list, tuple, set !")
        if not len(vec) == 3:
            raise IndexError("Your vector should have len = 3 !")
        vec_nump = numpy.array([vec])
        assert vec_nump.shape == (1,3), \
                "Wrong input: vec should be a nonnested container of lengtn 3!"
        return vec_nump
    @classmethod
    def tests(Cls):
        #tests initialization
        print("\nTesting class {}...\n".format(Cls))
        print("Initialization...")
        obj = Cls(lattice = (1., 1., 1., 90., 90., 90.))
        print("PASSED !")
        try:
            obj1 = Cls(lattice = [1,2,3, 90, 90, 90])
        except TypeError:
            print("Wrong type of argument upon initialization PASSED !")
        try:
            obj2 = Cls(lattice = (1,2,3))
        except IndexError:
            print("Wrong number of arguments upon initialization: PASSED !")

        #print("Here is the lattice: "),obj.lattice
        print("Here is the tensor"), obj.tensor
        
        print("Testing distance method...")
        print("Distance: {}".format(obj.distance([1,0,0],[0,1,0])))
        print("PASSED !")

        print("Testing length method...")
        print("Length: {}".format(obj.length([1,0,0])))
        print("PASSED !")

        try:
            obj.length(numpy.array([[1,2,3]]))
        except TypeError:
            print("Wrong vector type test PASSED !")

        try:
            obj.length([1,2])
        except IndexError:
            print("Wrong number of vector components test PASSED !")

        try: 
            obj.length([[1,2],[3,4],[5,6]])
        except AssertionError:
            print("Wrong vector shape assertion test PASSED !")

        print("Testing inversion...")
        print("Here is the inverted object {}".format(obj.inverted().tensor))
        print("DONE WITH TESTING!")

def symopSymList(sgNr):
    """For a S.G. Returns list of symmetry operations in symbolic form"""
    sg_hall = sgtbx.space_group_symbols(sgNr).hall()
    sg = sgtbx.space_group(sg_hall)
    print "\nGENERATING LIST OF SYMMETRY OPERATIONS FOR S.G. %d..." % (sgNr)
    symop_list = []
    sg_order = sg.order_z()
    for i in range(0, sg_order):
        symop_list.append(sg(i).as_xyz())
    print "A LIST OF %d SYMMETRY OPERATIONS WAS CREATED!\n" % (sg_order)
    return symop_list

def symopMatList(sgNr):
    """Returns list of symmetry operations in matrix(sympy) form """
    return map(xyzt2augmat, symopSymList(sgNr))

def xyzt2augmat(xyzt):
    """receives parsed symmetry operations from the Parser and returns the
    corresponding augmented matrix"""
    if not isinstance(xyzt, str):
        return xyzt

    row_map = {
                    "x":[1,0,0,0], #"-x": [-1,0,0,0],
                    "y":[0,1,0,0], #"-y": [0,-1,0,0],
                    "z":[0,0,1,0]} #"-z": [0,0,-1,0]}

    aug_m = []
    shifts = []
    xyzt_list = Parser(xyzt)

    for item in xyzt_list:
        mult_str, letter, transl_str = item
        row_in = row_map[letter]
             
        if '/' in mult_str:
            numer_mult, denom_mult = mult_str.split('/')
            mult = sympy.Rational(eval(numer_mult),eval(denom_mult))
        else:
            mult = eval(mult_str)

        if '/'in transl_str:
            numer_transl, denom_transl = transl_str.split('/')
            transl = sympy.Rational(eval(numer_transl), eval(denom_transl))
        else:
            transl = eval(transl_str)
        row_out = [nr*mult for nr in row_in]
        row_out[-1] = transl
        aug_m.append(row_out)
    aug_m.append([0,0,0,1])
    return sympy.Matrix(aug_m)
        
def xyz2mat(xyz):
    """Given transformation in symbolic form returns corresponding sympy 
    matrix"""
    mydict = {"x":[1,0,0], "-x": [-1,0,0], 
              "y":[0,1,0], "-y": [0,-1,0],
              "z":[0,0,1], "-z": [0,0,-1]}
    x = [mydict[item] for item in [item.strip() for item in xyz.split(",")]]
    return sympy.Matrix(x)

def symop2numpy(sym_op):
    """Returns symop numpy array. symop should be either in
    symbolic form or an instance of sympy.Matrix. If it was in symbolic form 
    it will be first transformed into sympy.Matrix"""
    if isinstance(sym_op, str):
        augm_mx = xyzt2augmat(sym_op)
    elif isinstance(sym_op, sympy.Matrix):
        augm_mx = sym_op
    else:
        raise TypeError("Wrong type of symmetry operation!")
    #transforming sympy matrix into numpy array
    augm_arr = numpy.array(numpy.array(augm_mx), numpy.float)
    return augm_arr


def modulo_n(mat,n = 1):
    """Performs operation modulo n under the translation part of the 
    augmented mat, where shift vector is in the last column"""
    d_x = xyzt2augmat('x+1,y,z')
    d_y = xyzt2augmat('x,y+1,z')
    d_x_inv = xyzt2augmat('x-1,y,z')
    d_y_inv = xyzt2augmat('x,y-1,z')

    while mat[0,-1] >= n:
        mat = d_x_inv*mat

    while mat[0,-1] < 0:
        mat = d_x*mat

    while mat[1,-1] >= n:
        mat=d_y_inv*mat

    while mat[1,-1] < 0:
        mat=d_y*mat
    
    return mat

def augmat2xyzt(augmat):
    """Returns symbolic representation of the augmented matrix"""
    index_symb = {0:'x',1:'y',2:'z'}
    columns = {
                "x" : augmat[0,0:3], 
                "y" : augmat[1,0:3],
                "z" : augmat[2,0:3],
                "shift" : augmat[0:3,3]
                }
    columns_parsed = {
                        "x":"", 
                        "y":"",
                        "z":"",
                        "shift":["","",""]
                        }

    for key, column in columns.iteritems():
        if key != "shift":
            for index,value in enumerate(column):
                if value > 0:
                    if columns_parsed[key] == "":
                        columns_parsed[key] += str(index_symb[index])
                    else:
                        columns_parsed[key] += "+" + str(index_symb[index])
                elif value < 0:
                    columns_parsed[key] += "-" + str(index_symb[index])
        else:
            for index, value in enumerate(column):
                if value > 0:
                    columns_parsed[key][index] ="+" + str(value)
                elif value < 0:
                    columns_parsed[key][index] = str(value)

    x = columns_parsed["x"] + columns_parsed["shift"][0] + ","
    y = columns_parsed["y"] + columns_parsed["shift"][1] + ","
    z = columns_parsed["z"] + columns_parsed["shift"][2] 

    string_rep = x + y + z
    return string_rep

    

def cosets_mod(stab_fragment, stab, n = 1):
    """Performs left coset decomposition stab:stab_fragment; 
    n-translation modulo. Returns list of lists[[coset1], [coset2], ...]"""
    cosets = []
    used_operations = []
    for ind,symop in enumerate(stab):
        if symop not in used_operations:
            l_coset = [modulo_n(symop*item, n) for item in stab_fragment]
            cosets.append(l_coset)
            used_operations.extend(l_coset)
    return cosets
                                            
def eqv_frag_ind(f1,f2,l1,l2):
    """Receives two cosets of the fragments and cosets of the layers
    in which the fragments are and returns a tuple of integers with the
    numbers designating position of fragments in those layers"""

    l1_bool = []
    l2_bool = []
    for i, cset in enumerate(l1):
        bool_array = [item in cset for item in f1]
        l1_bool.append(all(bool_array))
    
    if any(l1_bool):
        for j, cset in enumerate(l2):
            bool_array = [item in cset for item in f2]
            l2_bool.append(all(bool_array))

    else:
        l1_bool = []
        l2_bool = []

        for i, cset in enumerate(l1):
            bool_array = [item in cset for item in f2]
            l1_bool.append(all(bool_array))

        for j, cset in enumerate(l2):
            bool_array = [item in cset for item in f1]
            l2_bool.append(all(bool_array))

    try:
        i_ind = l1_bool.index(True)
        j_ind = l2_bool.index(True)

        return i_ind, j_ind
    except ValueError:
        print "Here are the fragments"
        print f1
        print f2
        raise ValueError

def cosets_are_disjoint(coset_list):
    """Checks that cosets from the coset_list are disjoint
    Probably quite slow step, because sympy.Matrix which are in cosets are
    unhashable, so can not use set-type"""


    symops = []
    for coset in coset_list:
        symops.extend(coset) 
    counts = []
    for item in symops:
        counter = -1
        for another in symops:
            if item == another:
                counter += 1
        counts.append(counter)
    if any(counts):
        return False
    return True

def is_group(symop_list, modulo = 1):
    """Checks whether symop_list constitute a group modulo in xy direction.
    IMPORTANT: only closure, existence of inverses and unit element are
    are checked"""
    #1.Closure

    #2.Inverses

    #3.Unit el-t
    pass

def shiftOrigin(symop = 'x,y,z'):
    """Shifts the origin of the group with symmetry operations in symopsList"""
    T_o =  xyzt2augmat(symop) #origin shift
    T_o_inv = T_o.inv()
    return lambda symopsList: map(lambda x: T_o*x*T_o_inv, 
                                  map(xyzt2augmat, symopsList))

def flatMap(inpList):
    """Flattens the list of lists"""
    outList = []
    for item in inpList:
        outList.extend(item)
    return outList


def applySymops(shifts = None):
    """Applyes shifts to symops from symopsList"""
    if shifts:
        return lambda symopsList: flatMap(
                    map(lambda shift: map(
                            lambda symop: shift*symop, map(xyzt2augmat,shifts)),
                        map(xyzt2augmat,symopsList)))
    

#left coset multiplication
lcoset = lambda (alist, symop): [symop*item for item in alist] 
#right coset multiplication
rcoset = lambda (alist, symop): [item*symop for item in alist]
#stabilizer after left coset multiplication
lstab = lambda (alist, symop): [symop*item*symop.inv() for item in alist]
#checks that coset is in stabilizer
coset_in_stab = lambda (coset, stab): all([symop in stab for symop in coset])
#checks that stabilizer of a subgroup devides that of a group
divides = lambda(H, G): len(G) % len(H) == 0


if __name__ == "__main__":
    MetricTensor.tests()


    shifts_x = map(xyzt2augmat, ['x,y,z', 'x+1/3,y,z','x+2/3,y,z'])
    shifts_y = map(xyzt2augmat, ['x,y,z', 'x,y+1/3,z','x,y+2/3,z'])
    coset_reps = map(xyzt2augmat, ['x,y,z', '-x,y,z', 'x,-y,z', '-x,-y,z'])

    shifts = []
    symops = []
    
    cosets_0 = []
    for coset in coset_reps:
        for s_x in shifts_x:
            for s_y in shifts_y:
                tr = s_y*s_x
                cosets_0.append(modulo_n(tr*coset, 1))

    shift_0 = sympy.Matrix([[0],[0],[0],[0]])
    shift_x_1_3 = sympy.Matrix([[sympy.Rational(1,3)],[0],[0],[0]])
    shift_x_2_3 = sympy.Matrix([[sympy.Rational(2,3)],[0],[0],[0]])
    
    # generating all shifts
    for shy in shifts_y:
        shifts.extend(map(modulo_n, lcoset((shifts_x, shy))))

    # generating all coset representative for the zeroth layer
    for shift in shifts:
        symops.extend(map(modulo_n, lcoset((coset_reps, shift))))

