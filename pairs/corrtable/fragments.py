import copy
import sys

from pairs.symmlib import core

class Fragment(object):
    def __init__(self, stab = [], coset = [], origin_shift = 'x,y,z'):
        """Fragment instance can be initialize in 2 ways:

        1.Fragment(stab)  --with optional origin_shift
            In this case it is assumed that given fragment is basis, i.e
            it will be used for generating further fragments by applying sym.
            operations. In this case:
            fragment.coset =  fragment.stab
        2.Fragment(stab, coset) --with optional origin_shift
            In this case it is assumed that this fragment was derived from
            another one by applying sym. operation, say g, from the left 
            In this case both coset and stab have to be provided. They are
            computed as follows
            stab = g*other.stab*g.inv()
            coset =  g*other.coset,
            where other- is the fragment from which the current one is derived
        """
        self._stab = []
        self._coset = []

        self.__initialize(stab, origin_shift, coset)

    def __initialize(self, stab, origin_shift, coset):
        """Initialization of the instance of with the user provided 
        parameters """
        self.stab = stab
        self.shift_stab_origin(origin_shift)
        if not coset:
            self.coset = copy.deepcopy(self.stab)
        else:
            self.coset = coset
###############################################################################
    @property
    def stab(self):
        return self._stab
    @stab.setter
    def stab(self, stab):
        """See doc string for __init__"""
        stab_list = []
        for item in stab:
            if str == type(item):
                stab_list.append(core.xyzt2augmat(item))
            else:
                stab_list.append(item)
        self._stab = stab_list
###############################################################################
    @property
    def coset(self):
        return self._coset
    @coset.setter
    def coset(self, coset = []):
        coset_list = []
        for symop in coset:
            if str == type(symop):
                coset_list.append(core.xyzt2augmat(symop))
            else:
                coset_list.append(symop)
        self._coset = coset_list

    def is_basic(self):
        """Basic fragment is the one whose coset contains 'x,y,z' operation
        for it is used to generate other fragments"""
        return core.xyzt2augmat("x,y,z") in self.coset
        

    def shift_stab_origin(self, origin_shift):
        """Shifts origin of the stabilizer of the instance
        IMPORTANT!!! origin_shift is not a modulo operation!!!"""
        if str == type(origin_shift):
            originShift = core.xyzt2augmat(origin_shift)
        else:
            originShift = origin_shift

        items = [] 
        for item in self.stab:
            items.append(originShift*item*originShift.inv())
        self.stab = items


    @classmethod
    def leftmult(Cls, self, symop):
        """Returns instance of Fragment obtained from the current instance by
        left multiplication with symop. The initialization of the new
        instance is done according to case 2 described in __init__"""
        if str == type(symop):
            sym_mat = core.xyzt2augmat(symop)
        else:
            sym_mat = symop 
        
        #I am not sure how the triple product modulo should be executed!!!!!
        
        #new_stab = map(
        #        core.modulo_n, core.lstab((self.stab, sym_mat))
        #        )
        new_stab = []
        sym_mat_inv = core.modulo_n(sym_mat.inv(), 1)

        for item in self.stab:
            new_elt = core.modulo_n(core.modulo_n(sym_mat*item,1)*sym_mat_inv,1)
            new_stab.append(new_elt)
        
        #new_stab = [core.modulo_n(core.modulo_n(sym_mat*item,1)*sym_mat.inv(),1) for item in  self.stab]


                
        new_coset = map(core.modulo_n, core.lcoset((self.coset, sym_mat)))

        return Cls(new_stab, new_coset)

    @classmethod
    def self_test(Cls):
        """Some basic testing"""
        basic_fragment = Cls(**{
                           "stab": ['x,y,z', '-y,-x,-z'],
                           "origin_shift": 'x,y,z+1/4'
                            })

        fragment1 = Cls.leftmult(
            basic_fragment,
            '-x,-y,z')
        print basic_fragment.stab
        print
        print basic_fragment.coset
        print
        print fragment1.stab
        print
        print fragment1.coset
    
class Layer(Fragment,object):
    def __init__(self, stab = [], coset = [], fragments = [],
                    origin_shift = "x,y,z"):
        """Layer instance can be initialized in 2 ways:

        1.Layer(stab, fragments)  --with optional origin_shift
            In this case it is assumed that given Layer is basic, i.e
            it will be used for generating other layers by applying sym.
            operations. In this case:
            self.coset =  self.stab
            >If fragments contains only single fragment and 
                |fragments[0].stab|<=|self.stab|, it is assumed
                that Layer should generate self.fragments by left
                coset decomposition of self.stab:fragments[0].stab
                >>Checks to be performed before left coset decomposition: 
                    a.|fragments[0].stab| devides |self.stab|
                    b.fragments[0].stab in self.stab
                >>Checks to be performed after left coset decomposition:
                    c.each of cosets self.fragments[i].coset in self.stab
                    d.|self.fragments[i].coset| divides |self.stabilizer|
                    e.cosets: self.fragments[i].coset are disjoint
                    f.|self.fragments[i].stab| divides |self.stab|
                    g.self.fragments[i].stab in self.stab


            >If fragment contains multiple fragments it is assumed that:
                self.fragments = fragments
                >>Checks to be performed(*):
                    e.cosets: self.fragments[i].coset are disjoint
                    f.|self.fragments[i].stab| divides |self.stab|
                    g.self.fragments[i].stab in self.stab

                                                                

        2.Layer(stab, coset, fragments) --with optional origin_shift
            In this case it is assumed that given layer was generated from
            the other one by left multiplication by say g:
            self = Layer.leftumlt(other, g).
            So:
            self.fragments = fragments
            self.stab = stab
            self.coset = coset
            where stab =  g*other.stab*g.inv()
                  coset = g*other.coset
                  fragment[i] = Fragment(other.fragment[i], g)
            >>Checks to be performed:
                ...the same as in (*)
                        """
        
        super(Layer, self).__init__(stab, coset, origin_shift)
        self._fragments = []
        self._basicFragment = []
        self.__initialize(fragments)
        self.__shift_classes = None

    def __initialize(self, fragments):
        """Initializes instance from the corresponding
        user-provided values"""
        self.fragments = fragments
###############################################################################
    @property
    def fragments(self):
        return self._fragments
    @fragments.setter
    def fragments(self, fragments):
        nr_frags = len(fragments)
        assert nr_frags >= 0
        "The layer contains no fragments!"#fragments is non-empty
        
        if nr_frags == 1:
            item_fragments = []
            basic_fragment = fragments[0]
            order_basic = len(basic_fragment.stab) 
            order_self = len(self.stab)
            #basic checks that can be switched off
            #1a in __init__
            assert order_self % order_basic == 0,\
                "Order of Fragment.stab does not divide that of Layer.stab!"
            #1b ...
            assert all([item in self.stab for item in basic_fragment.stab]),\
                "There are symops in Fragment.stab that are not in Layer.stab!"


            lcosets = core.cosets_mod(basic_fragment.stab, self.stab)
            for lcoset in lcosets:
                item_fragments.append(Fragment.leftmult(basic_fragment, 
                                                                lcoset[0]))
            self._fragments = item_fragments
            #basic checks to do after coset decomposition
            #1c
            assert all([core.coset_in_stab((coset, self.stab)) \
                    for coset in lcosets]),\
                """Cosets of the consituent fragments
                        do not belong to the stabilizer of the layer!"""
            #1d
            assert all([core.divides((coset, self.stab)) for coset in lcosets]),\
                """Order of cosets of consituent fragments does not devides
                        order of the layers's stabilizer"""
        else:
            self._fragments = fragments
        #basic checks that can be switched off    
        coset_list = [fragment.coset for fragment in self._fragments]
        stab_list = [fragment.stab for fragment in self._fragments]
        #1e
        assert core.cosets_are_disjoint(coset_list), "Cosets are not disjoint!"
        #1f
        assert all([core.divides((stab, self.stab)) for stab in stab_list]),\
                """Order of stabilizer of some constituent fragments
                        does not devide order of the layers's stabilizer"""
        #1g
        assert all([core.coset_in_stab((stab, self.stab)) \
                    for stab in stab_list]),\
                """Stabilizers of the consituent fragments
                        do not belong to the stabilizer of the layer!"""
        
###############################################################################
    @property
    def basicFragment(self):
        """Returns instance of the basic fragment of self"""
        if self._basicFragment:
            return self._basicFragment

        for fragment in self.fragments:
            if fragment.is_basic():
                self._basicFragment = fragment
        return self._basicFragment
###############################################################################
    def sort_fragments(self, coset_reps):
        """Sorts list self.fragments. Current algorithm is quite primitive
        and requires list of cost_reps. The list of self.fragments is sorted
        based on the order of symmetry elements in coset_reps, which are
        coset representatives of the decomposition 
        self.stab:fragment.stab, where fragment.stab has to be the fragment
        which was used for generating the whole layer"""

        indices = ["*" for i in range(len(self.fragments))]
        #in which coset is symio
        for i,symop in enumerate(coset_reps):
            for j,fragment in enumerate(self.fragments):
                if symop in fragment.coset:
                    indices[i] = j
        fragments_s = []
        for ind in indices:
            fragments_s.append(self.fragments[ind])

        self.fragments = fragments_s

    def fragment_index(self, coset):
        """Returns index of the fragment that has coset attribute or None
        if there is no such a fragment."""
        set_length = len(self.fragments[0].coset)
        assert len(coset) == set_length,\
                "Wrong number of sym.ops provided. It has to be %d symops!" \
                %(set_length)
        for nr, fragment in enumerate(self.fragments):
            if all([symop in fragment.coset for  symop in coset]):
                return nr
        return None              
        
    def build_shift_classes(self):
        """ IMPORTANT: THIS METHOD IS PARTLY HARDCODED, FOR IT ONLY
        ACCOUNTS FOR SHIFTS WITHIN XY PLANE REQUIRED FOR MODELLING
        DISORDER IN ZEOLITE-BETA. IT CAN HOWEVER BE EASILY MODIFIED!!!
        --------------------------------------------------------------
        This method will group fragments into shift classes(which will
        be stored in self.__shift_classes)  in the
        following way:
            1. Find basic fragment within the layer
            2. Find fragments that are obtained from the basic one without
                applying shifts and put them into
                one class, shift_class_0. This is achieved by sear-
                ching for symops in self.stab that have zero shift component,
                applying them to basic_fragment.stab and locating indices of 
                the corresponding fragements within self.
            3.Shift_class_i is produced by applying one of all the possible
              shift vectors (symops of self.stab that have unit rotation
              matrix as their rotation part) to the fragments in 
              zero_group_0, locating so obtained fragments within the layer
              and storing them within as a single group.
            4.The total number of shift_classes so obtained will be equal
              to the number of the possible pure shifts.
            5.Those shift_classes are later employed for the computation of 
              shift vector between two fragments 
            """
            
        self.__shift_classes = []
        shift_class_0 = set()
        basic_fragment = self.basicFragment
        basic_stab = basic_fragment.stab
        all_shifts = set()

        unitAugMat = core.xyzt2augmat('x,y,z')
        zero_shift_xy = unitAugMat[0:2,-1] 
        unit_rot_part = unitAugMat[0:3,0:3]
        
        #Constructing shift_class_0 and filling all_shfits with indices
        #of pure shift symops in self.symop 
        for nr,symop in enumerate(self.stab):
            if symop[0:2,-1] == zero_shift_xy:
                coset = map(core.modulo_n, core.lcoset((basic_stab, symop)))
                index = self.fragment_index(coset)
                shift_class_0.add(index)
            elif symop[0:3, 0:3] == unit_rot_part:
                all_shifts.add(nr)

        self.__shift_classes.append(shift_class_0)
        

    @classmethod
    def leftmult(Cls, self, symop):
        """to be filled"""
        if str == type(symop):
            sym_mat = core.xyzt2augmat(symop)
        else:
            sym_mat = symop 
        temp = super(Layer,Cls).leftmult(self, symop)
        
        fragments = []

        if self.fragments:
            for fragment in self.fragments:
                fragments.append(Fragment.leftmult(fragment, sym_mat))

        kwargs = {
                    "stab": temp.stab,
                    "coset": temp.coset,
                    "fragments": fragments
                }
                                
        return Cls(**kwargs)

