import copy
import sys

from pairs.symmlib import core

class Bilayer(object):
    def __init__(self, layers = [],):
        """The instance of bilayer has to initialized in one way only:
            Bilayer(layers)
            len(layers) == 2
            one of the layers in the bilayer has to be generated from the other
            by some symmetry operation g; this is important for the correct
            calculation of the stabilizer ot the bilayer! (see self.stab)
            So the run checks are the following:
            a.layers[0].coset and layers[1].coset are disjoint
            b.|layer[0].coset| == |layer[1].coset|"""
        assert len(layers) == 2, "Wrong number of layers. Only 2 are allowed!"
        self.layers = layers
        self._stab = []

    @property
    def stab(self):
        """Computing stabilizer of a bilayer(or any 2 fragments).
        self.layers must contain 2 Layers L_0 and L_1. 

        IMPORTANT: L_1 = g*L_0 so that if  S_0 = stab(L_0) => L_1 = g*S_0*L_0. 
        The latter means tha coset g*S_O contains all symops to generate L_1
        from L_0.

        The stabilizer of the bilayer is then a union of 2 sets:
            1. set1 contains all common symmetry elements of S_0 and S_1. 
                El-ts of this set simply transform each of the layers into
                itself
            2. set2 contains symmetry elements that transform L_0 into L_1
               and at the same time L_1 into L_0. El-ts that transform L_0
               into L_1 are in the coset:
                g*S_0       (that's why it is important
                            that L_1 = g*L_O, that the corresponding coset
                            attributes of the layers are set correctly)
               so g' from g*S_0 applied to L_1 should yield L_0:
                g'*L_1=L_0 and since g'L_0 = L_1 =>g'*g'L_0 = L_0 => g'*g'
                or g'^2 should belong to S_0.
                Hence set2 contains elements of coset g*S_0 (or L_1.coset)
                which squares belong to S_0.
                """
        if self._stab:
            return self._stab
        #finding overlap between stabilizers of the constituent layers
        stab_overlap = \
                [item for item in self.layers[1].stab
                        if item in self.layers[0].stab]
        #for a layer that is generated from the other one by a stacking oper:
        #finding coset members whose 2 orders belong to stabilizer of the basic
        #layer
        unit_element = core.xyzt2augmat('x,y,z')
        sec_ord_elts = []
        for nr,layer in enumerate(self.layers):
            if unit_element not in layer.coset:
                for symop in layer.coset:
                    if core.modulo_n(symop*symop,1) in self.layers[1-nr].stab \
                            and symop not in stab_overlap:
                            sec_ord_elts.append(symop)

        self._stab = stab_overlap + sec_ord_elts
        #performing some basic checks
        coset_union = self.layers[0].coset + self.layers[1].coset
        #print len(self.layers[0].coset)
        assert core.cosets_are_disjoint([self.layers[0].coset,
            self.layers[1].coset]), "Cosets are not disjoint!"
        
        assert len(self.layers[1].coset) == len(self.layers[0].coset),\
          "Coset attributes for the constituent layers are of different order!"
        return self._stab

    def indexed_fragments(self, index_tuple):
        """Returns reference to the fragments indexed by index_tuple:
            index_tuple[0] corresponds to the index of the fragment in
            the self.layers[0].fragments and index_tuple[1] to
            self.layers[1].fragments"""
        return (self.layers[0].fragments[index_tuple[0]],
                    self.layers[1].fragments[index_tuple[1]])


    @staticmethod
    def check_stab(stab):
        """This method performes check of the stabilizer of a bilayer.
        In its present version it only checks that two layers of beta,
        one basic and the one generateb by applying m_z to the basic layer
        have the stabilizer pmmm, including corresponding shiftst"""

        pmmm = core.symopMatList(47)
        d_x = map(core.xyzt2augmat,['x,y,z','x+1/3,y,z','x+2/3,y,z'])
        d_y = map(core.xyzt2augmat,['x,y,z','x,y+1/3,z','x,y+2/3,z'])
        
        pmmm_dx = []
        for symop in d_x:
            pmmm_dx.extend(map(core.modulo_n,core.lcoset((pmmm, symop))))
        pmmm_dx_dy = []
        for symop in d_y:
            pmmm_dx_dy.extend(map(core.modulo_n, core.lcoset((pmmm_dx, symop))))

        print "\n----------Checking stabilizer of the bilayer!--------------\n"
        print "The length of the generated stabilizer is: %d" %( len(pmmm_dx_dy))
        print "The length of the tested stabilizer is: %d" % (len(stab))
        
        counter1 = 0
        counter2 = 0
        for item in stab:
            if item in pmmm_dx_dy:
                counter1 += 1
            else:
                print "This item is not in the generated stabilizer!"
                print item
                print
        for item in pmmm_dx_dy:
            if item in stab:
                counter2 += 1
            else:
                print "This item is not in the tested stabilizer!"
                print item
                print
        print "\nRESULT:Generated and tested stabilizers are equal:", \
                                                        counter1 == counter2
        print "\n--------------------DONE!----------------------------------\n"
def stab_shift(stab_list,*args):
    """Applies shifts listed in args to the symmetry elements in stab_list"""
    stab_shifted = []
    for symop in stab_list:
        for d in args:
            symop_d = xyzt2augmat(d)*symop
            stab_shifted.append(symop_d)
    return stab_shifted

    

class Pair(Bilayer):
    def __init__(self, indices, bilayer):
        """The pair is identified by the tuple of indices (i,j) of the 
        constituent fragments. Those indices correspond to the indices of
        the fragments in the bilayer.layer[0-1].fragments. We designate
        them X_i and X_j.
        Important: bilayer.layer[1] should be generated 
        from bilayer.layer[0], i.e. the latter must contain the basic fragment!

        IMOPRTANT:To allow correct computation of the self.stab we will have to
        reset coset attributes for X_i and X_j to make them 
        X_i.coset = X_i.stab =S_i
        and X_j.coset = g*S_i (See Bilayer.stab to understand why it is
        important)
        so we deepcopy fragments so that we do not change coset attributes 
        of those of them that consitute the bilayer!
        """
        self.__init_fragments =  [bilayer.layers[0].fragments[indices[0]],
                bilayer.layers[1].fragments[indices[1]]]

        fragments_copy = map(copy.deepcopy, self.__init_fragments)
         
        self._reset_cosets(fragments_copy)
        
        super(Pair, self).__init__(fragments_copy)
        self._bilayer = bilayer
        self._cosets = []
        self._indEqvPairs = []

    def _reset_cosets(self, fragments):
        """For finding the stabilizer of the pair based on the algorithm 
        described in Bilayer.stab(it is inherited by the current class):
        The theory fo finding the stabilizer of the pair:
            The pair consists of fragments X_i and X_j, which coset attributes
            (i.e. X_i.coset),g_i*S_0 and g_j*S_O are in the following 
            relation with the basic fragment:
               (1) g_i*S_0*X_0 = g_i'*X_0 = X_i, g_i' in g_i*S_0
               (2) g_j*S_0*X_0 = g_j'*X_0 = X_j, g_j' in g_j*S_0,
                    S_0 = stab(X_0)

            from (1): X_0 = g_i'.inv*X_0
            Putting that into (2): yields g_j'*g_i'.inv*X_i = X_j
            Hence: reset coset attributes to the following:            
                1. X_j.coset = g_j*X_i.coset.inv()
                2. X_i.coset = X_i.stab
                """
        XiCosInv = [core.modulo_n(symop.inv(), 1)\
                for symop in fragments[0].coset]
        g_j = fragments[1].coset[0]
        fragments[1].coset = [core.modulo_n(g_j*symop, 1)\
                for symop in XiCosInv]
        fragments[0].coset = fragments[0].stab
        #return fragments

    
    @property
    def cosets(self):
        """decomposing group of a bilayer with respect to the pair
        stabilizer. Coset representatives generate pairs, symmetry
        equivalent to self. Returns list of cosets of the decomposition"""
        if self._cosets:
            return self._cosets
        self._cosets = core.cosets_mod(self.stab, self.bilayer.stab, 1)
        return self._cosets

    @property
    def indOfEqvPairs(self):
        """Uses self.cosets to find indices of fragmens in bilayer that form
        pairs symmetry equivalent to the self. Indices are a a list of tuples:
        [(),...(j,i),...()], j-index of the fragment in self.bilayer.layers[1]
        and i-...in self.bilayer.layers[0]. Theory beyound:
            Coset representatives g of decomposition
            self.bilayer.stab:self.stab,
            stored in self.cosets produce pair symmetry equivalent to self.
            We designate self=X_iX_j. 
            (1) X_i = g_i*S_0*X_0
            (2) X_j = g_j*S_0*X_0
            Acting on (1) and (2) with g we get:
                g*X_i = g*g_i*S_0*X_0
                g*X_j = g*g_j*S_0*X_0
            So we only have to find indices of cosets g*g_i*S_0 and g*g_j*S_0
            in coset attributes of self.bilayer.layers[0/1]"""
        if self._indEqvPairs:
            return self._indEqvPairs
        #generating cosets of the equivalent pairs with respect to the
        #basic fragment stabilizer.
        indices = []
        for cset in self.cosets:
            gg_iS_0 = [core.modulo_n(cset[0]*symop, 1)\
                    for symop in self.__init_fragments[0].coset]
            gg_jS_0 = [core.modulo_n(cset[0]*symop, 1)\
                    for symop in self.__init_fragments[1].coset]
            #finding the indices of the cosets among the cosets of layer0
            # and layer1
            cosLay0Fragms = [fragment.coset\
                    for fragment in self.bilayer.layers[0].fragments]
            cosLay1Fragms = [fragment.coset\
                    for fragment in self.bilayer.layers[1].fragments]
            #finding indices of the corresponding cosets
            ind01_tuple = core.eqv_frag_ind(gg_iS_0, gg_jS_0,\
                                            cosLay0Fragms, cosLay1Fragms)

            indices.append(ind01_tuple)

        self._indEqvPairs = indices

        return self._indEqvPairs
    

    @property
    def bilayer(self):
        return self._bilayer

    def shift_vec(self):
        "Computes shift vector between 2 fragments of the pair"
        return self.layers[1].coset[0][:,-1] - self.layers[0].coset[0][:,-1]
        #return self.__init_fragments[1].coset[0][:,-1] -\
        #        self.__init_fragments[0].coset[0][:,-1]

