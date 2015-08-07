import copy
import sys

from pairs.symmlib import core
from pairs.corrtable.pairfragments import Pair

class PairCorrTab(object):
    """Pair correlation table like it is implementted in Yell.
    So for a bilayer containing layer0 and layer1 the table contains
    n x n entries, provides layer1 is produced from layer0 by some symop.
    and hence both layers have equivalent number of constituent fragments.
    """
    def __init__(self, bilayer):
        self._content = []
        self._bilayer =  bilayer
        self._nrUniqPairs = []
        self._uniqPairs = []
        self._symRelatedPairs = []
        self._corrValues = None

    @property
    def bilayer(self):
        return self._bilayer

    @property
    def content(self):
        """Complete pair correlation table where ich entry (j,i) corresponds
        to pair with indices (i,j), i-belonging to the first and j to the 
        second of bilayer.layers, each entry contains a number, so that symm.
        related pairs contain the same number,i.e if there are 26 symmetry
        independent pairs than this number will be from the set (1,2,...26)"""
        if self._content:
            return self._content
        self.corr_analysis()
        return self._content

    @property
    def uniqPairs(self):
        """Returns dict with {nr:(i,j),....)},
        for nr see self.content.__doc__. (i,j)is as single representative
        for a group of symmetry related pairs,i.e if there are 26 symmetry
        independent pairs than this number will be from the set (1,2,...26)"""
        if self._uniqPairs:
            return self._uniqPairs

        uniqPairs = dict()
        
        for k in self.symRelatedPairs.keys():
            uniqPairs[k] = self.symRelatedPairs[k][0]
        self._uniqPairs = uniqPairs
        return self._uniqPairs

    @property
    def symRelatedPairs(self):
        """Gives dict {nr:[(i,j),(i,j)],....}, where nr corresponds to the 
        entry (j,i) in self.content"""
        if self._symRelatedPairs:
            return self._symRelatedPairs
        self.corr_analysis()
        return self._symRelatedPairs

    @property
    def nrUniqPairs(self):
        if self._nrUniqPairs:
            return self._nrUniqPairs
        self.corr_analysis()
        return max(map(max, self.content))

    @property
    def corrValues(self):
        """Returns previously supplied pair correlation values"""
        return self._corrValues

    @corrValues.setter
    def corrValues(self, **corrValues):
        """Sets up pair correlation probabilities for certain pairs.
        The supplied argument should be a dictionary where keys are labels
        for certain pairs"""
        pass

    @property
    def pickedPairs(self):
        """See self.pick_pairs for the description"""
        return self._pickedPairs

    def pick_pairs(self, shift_vectors):
        """Shift vectors is a dictionary: {key: shift_vector}, where
        shift_vector--is a shift vector between two fragments that form a pair.
        The procedure walks through self.uniqPairs, checks their shift_vec
        attribute and if that is in shift_vector list, it puts the info of 
        this vector into a dict in the folowing way:
        {key: [...{nr:(i,j)}....]}, where key is from the supplied
        shift vectors, {nr:(i,j)} is described in docstring for uniqPairs
        The final result is assigned to self._pickedPairs
        """
        output = dict()
        for key1 in shift_vectors.keys():
            vecs = shift_vectors[key1]
            for key2 in self.uniqPairs.keys():
                pair = Pair(self.uniqPairs[key2], self.bilayer)
                shift = pair.shift_vec()[0:2]
                #print "YESS!"
                if shift in vecs:
                    output[key1] = \
                        output.get(key1,[])+[{key2:self.uniqPairs[key2]}]
        self._pickedPairs = output
    def subs(self, check_val, sub_val):
        """In self.content substitutes every entry with check_val for sub_val"""
        for j,line in enumerate(self.content):
            for i,entry in enumerate(line):
                if entry ==  check_val:
                    self.content[j][i] = sub_val
    def fill_with_values(self, kwargs):
        """kwargs = {entry_label: value_to_fill}. Fills table entries
        with the specified entry_labels with values_to_fill"""
        keys =  kwargs.keys()
        for nr1, line in enumerate(self.content):
            for nr2, entry in enumerate(line):
                if entry in keys:
                    self.content[nr1][nr2] = kwargs[entry]
                else:
                    self.content[nr1][nr2] = 0

    def corr_analysis(self):
        """Generating pair correlation table which i-th entry of the j-th line
        contains *, which later will be filled with instance of ji pair,
        j-th fragment coming from the zeroth and i-th fragment from some n-th
        layer of a bilayer"""
        
        pairCorrTable = []

        for ind1 in range(len(self.bilayer.layers[0].fragments)):
            tab_line=[]
            for ind2 in range(len(self.bilayer.layers[1].fragments)):
                tab_line.append("*")
            pairCorrTable.append(tab_line)
        counter = 1
        symRelatedPairs = dict()

        for j,line in enumerate(pairCorrTable):
            for i,entry in enumerate(line):
                #print i,j
                if entry == "*":
                    indices = Pair((i,j),self.bilayer).indOfEqvPairs
                    symRelatedPairs[counter] = indices # symeq pairs
                    for ind0,ind1 in indices:
                        if pairCorrTable[ind1][ind0] != "*":
                            print "Warning! The entry (%d,%d) is already in use" % (ind1, ind0)
                        try:
                            pairCorrTable[ind1][ind0] = counter
                        except IndexError:
                            print ind1,ind0
                            raise
                    counter += 1  
                    self._content = pairCorrTable

        self._nrUniqPairs = counter
                    
        self._symRelatedPairs = symRelatedPairs

        self._content = pairCorrTable

    def check_content(self):
        sums = dict()
        sum_rows = map(sum, self.content) 
        print sum_rows
        sum_columns = [0 for i in range(len(self.content[0]))]
        
        for row in self.content:
            for nr,item in enumerate(row):
                sum_columns[nr]+=item
        return sum_rows == sum_columns

    def check_corr(self):
        """This test should check that for every pair in the bilayer:
            symmetry operators of the stabilizer of a bilayer generate
            the same symmetry equivalent pairs as obtained by coset decom
            position of the stabilizer of the bilayer with respect to the sta
            bilizer of a pair(method Pair.indOfEqvPairs)"""
        symops = self.bilayer.stab
        frags_0 = self.bilayer.layers[0].fragments #fragments of the first layer
        frags_1 = self.bilayer.layers[1].fragments #fragments of the second layer

        cosets_0 = [getattr(fragment, "coset") for fragment in frags_0]
        cosets_1 = [getattr(fragment, "coset") for fragment in frags_1]
        alerts = []
        checks = []
        for i,fr0 in enumerate(frags_0):
            for j,fr1 in enumerate(frags_1):
                print "-----------Testing pair (%d,%d)-----------------------"\
                                                        % (i,j)

                #calculated indices for the symmetry equivalent pairs
                pairs_calc = Pair((i,j),self.bilayer).indOfEqvPairs
                
                
                #observed indices for the symmetry equivalent pairs
                pairs_obs = []
                for symop in symops:                 
                    cs0 = map(core.modulo_n, core.lcoset((fr0.coset,symop)))
                    cs1 = map(core.modulo_n, core.lcoset((fr1.coset,symop)))
                    indices = core.eqv_frag_ind(cs0,cs1,cosets_0, cosets_1)
                    pairs_obs.append(indices)
                print "Number of eqv. pairs calc.,obs.: %d, %d" \
                                            % (len(pairs_calc),len(pairs_obs))
                set_calc = set(pairs_calc) 
                set_obs = set(pairs_obs)
                print "Number of unique eqv. pairs calc.,obs.: %d, %d" \
                                % (len(set_calc),len(set_obs))
                
                test_value = all((set_obs <= set_calc, set_obs >= set_calc))
                checks.append(test_value)
                if not test_value:
                    alerts.append((i,j))
                print "Calculated and observed pairs are the same:", test_value
                sys.stdout.flush()
        print "\nALL THE PAIRS PASSED THE TEST:", all(checks)
        print "Problematic pairs (if any):\n"
        print alerts

    def __repr__(self):
        """Prints Pair Correlation Table"""
        output = ""
        nrOfLines = len(self.content)-1
        for nr,aline in enumerate(self.content):
            line = ",".join(map(str,aline))
            if nrOfLines != nr:
                output += line + ',\n'
            else:
                output += line 
        return output

