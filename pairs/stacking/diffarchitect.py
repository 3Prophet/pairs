import sys
import os
import re
import numpy
import copy
import sympy
import math

from iotbx.cif import reader
from iotbx.cif import model
from iotbx.shelx import writer
from iotbx.shelx import from_ins
from cctbx import xray
from cctbx import crystal
from cctbx import adptbx
from cctbx.array_family import flex
from cctbx import sgtbx

from pairs.symmlib import core

class Fragment(object):
    """This class is used to model a fragment of the structure"""
    def __init__(self, scatterers):
        """Initialization of a Fragment's instance with array of scatterers:
        arrray should be an instance of
        cctbx.array_family.flex.xray_scatterer(). Its members should
        be instances of cctbx.xray.scatterer """
        if not isinstance(scatterers, flex.xray_scatterer):
            raise TypeError("Initialize item with an instance of cctbx.array_family.flex.xray_scatterer")
        self._scatterers = scatterers

    def __iter__(self):
        return self._scatterers
    
    @property
    def scatterers(self):
        return self._scatterers

    def get_by_label(self, label):
        """Returns scatterere with the matching label"""
        LABEL = label.upper()
        if not isinstance(label, str):
            raise TypeError("Argument is not a string instance !")
        for scatterer in self._scatterers:
            if getattr(scatterer, "label").upper() == LABEL:
                return scatterer
        return None

    def get_index_by_label(self, label):
        """For the scatterer with a label, returns its
        index in self._scatterers"""
        LABEL = label.upper()
        if not isinstance(label, str):
            raise TypeError("Argument is not a string instance !")
        for index,scatterer in enumerate(self._scatterers):
            if getattr(scatterer, "label").upper() == LABEL:
                return index
        return None

    def gen_unique_labels(self):
        """Provides unique labels to each member of self.scatterers array"""
        sorted(self._scatterers, key = lambda x: x.label)
        counted = dict()
        for scatterer in self._scatterers:
            elt = scatterer.scattering_type
            counted[elt] = counted.setdefault(elt, 0) + 1
            setattr(scatterer, "label", "{0:s}{1:d}".format(
                    elt,
                    counted[elt])
                    )
    def limit_to_range(self, **kwargs):
        """Returns a copy of object with its scatterers falling into coordina
        te ranges specified by **kwargs. For example kwargs could be:
            x = (0,0.5), z = (0, 0.8)"""
        map_kwargs = {'x':0, 'y':1, 'z':2}
        init_scatterers = flex.xray_scatterer()
        for scatterer in self._scatterers:

            check_list =[[scatterer.site[map_kwargs[key]], kwargs[key]]  
                    for key in kwargs]
            cond = all([ all([item[0] >= item[1][0], item[0] <= item[1][1]]) 
                            for item in check_list]) 
            if cond:
                init_scatterers.append(scatterer)
        return Fragment(init_scatterers)
                            
    def set_value(self, **kwargs):
        """Sets an attribute for each member of self.scatterers
        to certain value, according to kwargs.  For a list of available 
        attribute please use help(cctbx.xray.scatterer)"""
        for scatterer in self._scatterers:
            for key,val in kwargs.items():
                setattr(scatterer, key, val)

    def set_occupancy(self, occupancy = 1, **kwargs):
        """Sets general multiplicity for scatterers. The multiplicity for the
        special scatterers are given in **kwargs,e.g. {'Si1':0.5}"""
        for scatterer in self._scatterers:
            scatterer.occupancy = kwargs.get(
                        scatterer.label,
                        occupancy
                        )

    def print_for_Yell(self, multiplicity = 1):
        """prints fragment in the format of Yell"""
        for scatterer in self._scatterers:
            print(
            "{0:4s} = {1:4s}{2:4f}{3:10.6f}{4:10.6f}{5:10.6f}{6:>10f}".format(
                            scatterer.label,
                            scatterer.scattering_type,
                            scatterer.occupancy,
                            scatterer.site[0],
                            scatterer.site[1],
                            scatterer.site[2],
                            scatterer.u_iso)
                )
    
    def deepcopy(self):
        """Returns deepcopy of a Fragment's instance"""
        return copy.deepcopy(self)

    def clone(self, sym_op):
        """Returns new instance of a Fragment transformed by sym_op"""
        new_frag = self.deepcopy()
        new_frag.transform(sym_op)
        return new_frag

    def transform(self, sym_op):
        """Mutates "site" attribute (coordinates) of each member of
        self.scatterers array with symm_op Example: symm_opp = 'x,y,z+1/4'"""
        #getting sympy matrix fo the symbolic sym_op
        augm_arr = core.symop2numpy(sym_op)

        for scatterer in self._scatterers:
            trans_xyz(scatterer, augm_arr)
            #print scatterer.site

    def transform_cond(self, sym_op = 'x,y,z', cond = 'x<0'):
        """applies sym_op to the scatterer whos x,y or z coordinate satisfies
        to condition cond"""
        xyz_map = {'x':0,'y':1,'z':2}
        parsed_cond = re.split(r"([<=>(<=)(>=)])",cond)
        ind = xyz_map[parsed_cond[0]]
        augm_mx = augm_arr(sym_op)

        for scatterer in self._scatterers:
            if eval(str(scatterer.site[ind])+''.join(parsed_cond[1:])):
                trans_xyz(scatterer, augm_mx)

    def set_value_atom(self, label, **kwargs):
        """Changes the value of the attribute of a particular scatterer
        designated with label"""
        for scatterer in self._scatterers:
            if scatterer.label == label:
                for attribute, value in kwargs.items():
                    setattr(scatterer, attribute, value)
    
    def set_attr_value(self, attr1, cond1, attr2, newvalue = '*1'):
        """Sets attribute2 to newvalue, based on cond1 for attribute 1.
        newvalue has to be string giving an opeartion or value.
        In the former case value of attr2 is changed by the operation the
        string provides attr1 has to be a tuple of strings in case the
        corresponding attribute is a tuple and its specific range is of
        interest, for example ("site", "[0,1]") or ("occupancy","")-for 
        a single parameter value """

        for scatterer in self._scatterers:
            cond = eval(str(getattr(scatterer, attr1[0]))+attr1[1]+cond1)
            #handles only the case of a single valued attribute 2
            if cond:
                if not set("*+-/").isdisjoint(newvalue) :
                    new_attr_value = eval(str(getattr(scatterer,
                                            attr2))+newvalue)
                else:
                    new_attr_value = eval(newvalue)
                setattr(scatterer, attr2, new_attr_value)
    
    def xyz_modulo(self, n = 1):
        """Applies modulo n operation to each of the coordinates"""
        for scatterer in self._scatterers:
            x,y,z = getattr(scatterer, "site")
            xyz_new = x%n, y%n, z%n
            setattr(scatterer, "site", xyz_new)

    @classmethod
    def merge(Cls, *others):
        """Returns new instance of the Fragment. The scatterers in the new
        fragment are obtained by merging two frag"""
        scatterers = flex.xray_scatterer()
        for other in others:
            for scatterer in other.scatterers:
                scatterers.append(scatterer)
        return Cls(scatterers)

def augm_arr(sym_op = 'x,y,z'):
    """returns augmented array corresponding to the symmetry operation sym_op"""
    mtr = sgtbx.rt_mx(sym_op)
    #constructing augmented matrix of the symmetry operation
    augm_arr = numpy.reshape(
                numpy.array(mtr.as_4x4_rational(), dtype = numpy.float),
                (4,4)
                )
    return augm_arr

def trans_xyz(scatterer, augm_arr):
    """Transforms coordinates of the scatterer according to the symmetry operat
    ion with the corresponding augmented array represented by augm_arr.
    scatterer should be the instance of  xray.scatterer class """
    #transforming sympy matrix into numpy array
    if not isinstance(augm_arr, numpy.ndarray):
        raise TypeError
    
    augm_mx =  augm_arr
    augm_vec = numpy.reshape(
                    numpy.array([scatterer.site[0],
                                scatterer.site[1],
                                scatterer.site[2],
                                1.0]),
                    (4,1)
                    )

    #performing the transformation
    augm_vec_new = numpy.dot(augm_mx, augm_vec)
    new_coords = tuple((
                    augm_vec_new[0,0],
                    augm_vec_new[1,0],
                    augm_vec_new[2,0])
                    )
            #assigning new coordinates to the scatterer
    setattr(scatterer, "site", new_coords)

class Model(Fragment):
    """Model, which is saved into the cif file """
    def __init__(self, symmetry = None, fragments = None, 
                    gen_unique_labels = False,
                    min_distance_sym_equiv = 0.0001):
        """Initialization of the Model's instance with symmetry (should be an
        istance of cctbx.crystal.symmetry) and list of fragments. Each fragment
        should be an istance of Fragment. min_distance_sym_equiv - if atoms are
        located close to the special position(within the given value) some of
        their coordinates will be rounded. Give this parameter some small value
        if you want to avoid this."""
        self._fragments = fragments
        self._symmetry = symmetry
        self._scatterers = flex.xray_scatterer()
        self._special_position_settings = crystal.special_position_settings(
                self._symmetry, min_distance_sym_equiv)
        if fragments:
            for fragment in fragments:
                for scatterer in fragment._scatterers:
                    self._scatterers.append(scatterer)

            if gen_unique_labels == True:
                self.gen_unique_labels()

        self._iucr_structure = None

    @property
    def fragments(self):
        return self._fragments
    @fragments.setter
    def fragments(self, value):
        raise NotImplementedError
    @fragments.deleter
    def fragments(self):
        raise NotImplementedError


    @property
    def symmetry(self):
        return self._symmetry
    @symmetry.setter
    def symmetry(self, value):
        raise NotImplementedError
    @symmetry.deleter
    def symmetry(self):
        raise NotImplementedError

    @property
    def scatterers(self):
        return self._scatterers
    @scatterers.setter
    def scatterers(self, value):
        raise NotImplementedError
    @scatterers.deleter
    def scatterers(self):
        raise NotImplementedError

    @property
    def iucr_structure(self):
        #self._fragments = [Fragment(self._scatterers)] #?
        self._iucr_structure = xray.structure(
                special_position_settings = self._special_position_settings,
                scatterers = self._scatterers
                )                
        return self._iucr_structure
    
    @iucr_structure.setter
    def iucr_structure(self, value):
        raise NotImplementedError
    @iucr_structure.deleter
    def iucr_structure(self):
        raise NotImplementedError

    def transform(self, symop, **kwargs):
        """Returns new model with new basis: b_new = symop*b_old
        and new coordinates: xyz_new = symop.inv()*xyz_old and the
        THE SAME SPACE GROUP ! WARNING! Take care: symmetry operations are not 
        transformed (at the moment) so this function in its current state 
        is to be used for expanding or contracting the unit cell!
        Symop should be either in symbolic form or be an augmented matrix
        of type sympy.Matrix"""
        if isinstance(symop, sympy.Matrix):
            if symop.shape != (4,4):
                raise IndexError("Wrong shape of the matrix. Correct: 4 x 4 !")
            symop_sympy = symop
        elif isinstance(symop, str):
            symop_sympy = core.xyzt2augmat(symop)
        else: 
            raise TypeError("Wrong type of the symmetry operation. See the docstring !")
        scatterers = copy.deepcopy(self._scatterers)
        fragment = Fragment(scatterers)
        fragment.transform(symop_sympy.inv())
        #print self._symmetry.unit_cell().parameters() #unit_cell().__dict__
        lattice = self._symmetry.unit_cell().parameters()
        pars = []
        pars.extend(lattice[0:3])
        angles = lattice[3:]
        pars.append(0)
        latt_old = sympy.Matrix(pars)
        latt_new = list(latt_old.T*symop_sympy)[0:3]
        latt_new.extend(angles)
        symmetry_new = copy.deepcopy(self._symmetry)
        setattr(symmetry_new, "unit_cell", tuple(latt_new))
        return self.__class__(symmetry_new, [fragment])

    def save_ins(self,fpath, wlength = 0.7100, sof = 0.003, nrcycles = 20 ):
        generator = writer.generator(self.iucr_structure, 
                        wavelength = 0.7100,
                        #temperature = 298,
                        overall_scale_factor = 0.003,
                        full_matrix_least_squares_cycles = 20
                        )
        print
        with open(fpath,'w') as fobject:
            for line in generator:
                fobject.write(line)
    
    def save_CIF(self, fpath):
        """Saving model into cif file"""
     
        cif_object = model.cif()
        cif_block = model.block()
        cif_object["BEA"] = cif_block
        space_group = sgtbx.space_group(self.iucr_structure.space_group())

        #unit cell:
        cell_pars = self.iucr_structure.unit_cell().parameters()

        cif_block["_cell_length_a"] = cell_pars[0]
        cif_block["_cell_length_b"] = cell_pars[1]
        cif_block["_cell_length_c"] = cell_pars[2]
        cif_block["_cell_angle_alpha"] = cell_pars[3]
        cif_block["_cell_angle_beta"] = cell_pars[4]
        cif_block["_cell_angle_gamma"] = cell_pars[5]

        space_group_type = self.iucr_structure.space_group_info().type()
        cif_block["_symmetry_cell_setting"] = \
                                    space_group.crystal_system().lower()
        cif_block["_symmetry_Int_Tables_number"] = space_group_type.number()
        cif_block["_symmetry_space_group_name_H-M"] = \
                                        space_group_type.lookup_symbol()
        #cif_block["_space_group.name_Hall"] = space_group_type.hall_symbol()

        symop_loop = model.loop(
            header = ("_symmetry_equiv_pos_as_xyz",)
                            )
        for symop_id, symop in enumerate(space_group):
            symop_loop.add_row(("'{}'".format(symop.as_xyz()),))

        struct_loop = model.loop(
            header = ("_atom_site_label",
                    "_atom_site_type_symbol",
                    "_atom_site_fract_x",
                    "_atom_site_fract_y",
                    "_atom_site_fract_z",
                    "_atom_site_U_iso_or_equiv"))

        for scatterer in self.iucr_structure.scatterers():
            struct_loop.add_row(
                        (scatterer.label,
                        scatterer.scattering_type,
                        scatterer.site[0],
                        scatterer.site[1],
                        scatterer.site[2],
                        scatterer.u_iso))
                        
    
        cif_block.add_loop(symop_loop)
        cif_block.add_loop(struct_loop)
        
        with open(fpath,'w') as fobj:
            fobj.write(cif_object.__str__())


class Initializer(object):
    """Initializes scatterers and symmetry either from the *.xyz or *.cif"""
    def __init__(self, fname = None, **kwargs):
        self.fname = fname
        self.scatterers = None
        self.symmetry = None
        self._initialize(**kwargs)

    @property
    def unit_cell(self):
        return self.symmetry._unit_cell

    def _initialize(self, **kwargs):
        base = os.path.basename(self.fname)
        if base.endswith(".xyz"):
            self._initialize_from_XYZ(**kwargs)
        elif base.endswith(".cif"):
            self._initialize_from_CIF(**kwargs)
        elif base.endswith(".ins") or base.endswith(".res"):
            self._initialize_from_ins(**kwargs)


    def _initialize_from_CIF(self, **kwargs):
        iucr_structure = read_CIF(self.fname, **kwargs)
        self.scatterers = iucr_structure.scatterers()
        self.symmetry = iucr_structure.crystal_symmetry()

    def _initialize_from_XYZ(self, **kwargs):
        coordinates_list = read_XYZ(self.fname, **kwargs)
        self.scatterers = flex.xray_scatterer()
        for line in coordinates_list:
            line_list = line.strip().split()
            self.scatterers.append(
                xray.scatterer(
                    label = line_list[1],
                    scattering_type = line_list[0],
                    site = (
                    float(line_list[2]),
                    float(line_list[3]),
                    float(line_list[4]),
                    ),
                    **kwargs
                    ))

    def _initialize_from_ins(self, **kwargs):
        iucr_structure = from_ins.from_ins(self.fname, **kwargs)
        self.scatterers = iucr_structure.scatterers()
        self.symmetry = iucr_structure.crystal_symmetry()
        for scatterer in self.scatterers:
            setattr(scatterer, "u_iso", 
                    adptbx.u_star_as_u_iso(self.unit_cell, scatterer.u_star))

    def generate_fragment(self):
        return Fragment(self.scatterers)

    def generate_model(self, **kwargs):
        return Model(symmetry = self.symmetry,
                fragments = [Fragment(self.scatterers)], **kwargs)        
    
    @classmethod
    def tests(Cls):
        initializer = Cls(os.path.join(os.getcwd()+r"/test_files", "beta.res"))
        print initializer.symmetry._unit_cell
        for scatterer in initializer.scatterers:
            print getattr(scatterer,"label"), getattr(scatterer,"site"),\
                    adptbx.u_star_as_u_iso(initializer.symmetry._unit_cell,
                            scatterer.u_star)
        

def read_XYZ(fname, **kwargs):
    """reads coordinates from the file and"""

    patt = re.compile(r"(Elmt\s+Label.+)")

    fobj = open(fname, "r")
    content = fobj.read()
    fobj.close()

    match = re.search(patt, content)
    coordinates = content[match.end():]
    coordinates_list = coordinates.strip().split("\n")
    return coordinates_list

def read_CIF(fname, **kwargs):
    """Returns model read from CIF file"""
    print("\nREADING STRUCTURAL INFORMATION FROM {} ...".format(
                os.path.basename(fname)))
    iucr_structure = \
            reader(fname).build_crystal_structures()["BEA"]
    if kwargs.get("expand_to_p1", 0): 
        print "\n...EXPANDING SRTUCTURE FROM .CIF IN P1...\n"
        iucr_structure = iucr_structure.expand_to_p1(sites_mod_positive = True)
    print("...DONE!\n")
    return iucr_structure

if __name__ == "__main__":
    Initializer.tests()

