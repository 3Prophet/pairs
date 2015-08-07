import diffarchitect
import os
import core
from cctbx import sgtbx
import threading
import copy
import fnmatch
import multiprocessing
import sympy
import numpy

class SenderAdapter(object):
    """Adapter class emulates sender interface for testing purposes"""
    class AdapterFile(object):
        """Implements get() for Adapters' attributes"""
        def __init__(self, fname):
            self.fname = fname
        def get(self):
            return self.fname

    def __init__(self, path_dir, Ins_filename):
        self.path_dir = self.__class__.AdapterFile(path_dir)
        self.Ins_filename = self.__class__.AdapterFile(Ins_filename)
        self._observers = []
        

    def subscribe(self, observer):
        self._observers.append(observer)

    def run(self):
        """To emulate refinement process:
        We open .ins file, change coordinates of 2 atoms
        and save .res file"""
        dir_path = self.path_dir.get()
        ins_file = self.Ins_filename.get()
        print dir_path, ins_file
        res_file = ins_file.rstrip(".ins") + ".res"
        initializer = diffarchitect.Initializer(os.path.join(dir_path, ins_file))
        model_original = initializer.generate_model(gen_unique_labels = False)
        model_refined = copy.deepcopy(model_original)
        scatterer1 = model_refined.scatterers[0]
        print("Refining xyz of scatterer {}".format(getattr(scatterer1, 'label')))
        x1,y1,z1 = getattr(scatterer1, 'site')
        y1 += 0.002
        setattr(scatterer1,'site', (x1,y1,z1))
        model_refined.scatterers[0] = scatterer1
        scatterer2 = model_refined.scatterers[1]
        print("Refining xyz of scatterer {}".format(getattr(scatterer2, 'label')))
        x2,y2,z2 = getattr(scatterer2, 'site')
        y2 += 0.003
        setattr(scatterer2,'site', (x2,y2,z2))
        model_refined.scatterers[1] = scatterer2
        model_refined.save_ins(os.path.join(dir_path, res_file))
        
        if self._observers:
            for observer in self._observers:
                observer.update(self)
        else:
            print("Message from: {}. No observers to update !".format(self.__class__.__name__))

    @classmethod
    def tests(Cls, testdir):
        dirpath = testdir
        for dir_path, dirnames, fnames in os.walk(testdir):
            dirnames = []
        inss = fnmatch.filter(fnames, "*.ins")
        if len(inss) == 0:
            raise IOError("No .ins in the directory {}".format(dir_path))
        if len(inss) > 1:
            raise IOError("More than one .ins file in the directory {}".format(dir_path))
        
        fname = inss.pop()
        sender = Cls(dirpath, fname)
        updater = ModelUpdater(os.path.join(testdir, "BEA_fragment_cut.cif"),
                "1/3x, 1/3y, 1/2z",
                1)
        sender.subscribe(updater)
        sender.run()
       
class RefinementObserver(object):
    """This class follows the refinement done by SHELX and updates
    its model_ins and model_res attributes based on the models contained
    in *.ins and *.res files after the refinement"""
    def __init__(self):
        self._model_ins = None #model contained in .ins file
        self._model_res = None #model contained in .res file
        self._observers = []

    @property
    def model_ins(self):
        return self._model_ins
    @model_ins.setter
    def model_ins(self, value):
        raise AttributeError
    @model_ins.deleter
    def model_ins(self):
        raise AttributeError

    @property
    def model_res(self):
        return self._model_res
    @model_res.setter
    def model_res(self, value):
        raise AttributeError
    @model_res.deleter
    def model_res(self):
        raise AttributeError

    def subscribe(self, observer):
        self._observers.append(observer)

    def _update_model(self, model_path, message):
        if os.path.exists(model_path) and os.path.isfile(model_path):
            initializer = diffarchitect.Initializer(model_path,  
                    min_distance_sym_equiv = 0.0001)
            print(message)
            return initializer.generate_model(gen_unique_labels =  False)
        else:
            print("No filename {} found !!!". format(model_path))
            return None
                   
    def _update_model_ins(self, path_ins):
        self._model_ins = self._update_model(path_ins,
                            "UPDATING INITIAL MODEL DONE !")
       
    def _update_model_res(self, path_res):
        self._model_res = self._update_model(path_res,
                            "UPDATING REFINED MODEL DONE !")
 
    def update(self, sender):
        """In the present version sender is an instance of flipGUI class"""
        refinement_dir = sender.path_dir.get()
        fname_ins = sender.Ins_filename.get()

        if not fname_ins.endswith(".ins"):
            fname_ins += ".ins"   
        path_ins = os.path.join(refinement_dir, fname_ins)
        path_res = os.path.join(refinement_dir, fname_ins.rstrip(".ins")+".res")
        #self._update_model_ins(path_ins)
        #self._update_model_res(path_res)

        threads = [
                threading.Thread(target = self._update_model_ins,
                                    args = (path_ins,)),
                threading.Thread(target = self._update_model_res,
                                    args = (path_res,))
                ]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        if self._observers:
            for observer in self._observers:
                multiprocessing.Process(target = observer.update,
                                        args = (self,)).start()

class ModelUpdater(RefinementObserver):
    def __init__(self, model_path, symop, modulo_n, instructions = None,
                                            **kwargs):
        """ Assume we have a model that we have to transform somehow, refine
        and than transform back the refined model, to get the improved version
        of the original model. This class does the job.
        Upon initialization it accepts: 
        >model_path -- an instance of the class Model
        >symop -- symmetry operation to transform initial model to the model
        that will be refined
        The self.model_updated will contain updated(after refinement) version 
        of the initial model.
        >instructions -- if during the refinement you change coordinates of
        your atoms manualy by means of some symmetry operation, than you should
        provide this symop together with the label of the atom so that
        the class can use it to get back the coordinate of the atom. This param
        is dict of type {..., atom_label: symop,...}, where symop is either
        symbolic or augmented sympy matrix.         
        """
        super(ModelUpdater, self).__init__()

        if not isinstance(model_path, str):
            raise TypeError("Path to the model file is not an instance of str !")
        if os.path.exists(model_path) and os.path.isfile(model_path):
            initializer = diffarchitect.Initializer(model_path)
            self._model = initializer.generate_model(**kwargs)
        if not isinstance(modulo_n, int):
            raise TypeError("Modulo operator is not an instance of int !")
        self._modulo_n = modulo_n
        self._model_path = model_path
        self._instructions = instructions 

        if isinstance(symop, sympy.Matrix):
            self._symop = symop
        elif isinstance(symop, str):
            self._symop = core.xyzt2augmat(symop)
        else:
            raise TypeError("Symm.op. is not an instance of sympy.Matrix or str !")
        self._model_intermediate = self._create_model_intermediate()
        self._model_updated = copy.deepcopy(self._model)

    def _create_model_intermediate(self):
        """The sole purpose of this function is to be able to maintain
        the whole part of each of the coordinates chopped by modulo 1
        operation."""
        return self._model.transform(self._symop)

    def _update_model_updated(self):
        for scatterer in self.model_res.scatterers:
            x0, y0, z0 = getattr(scatterer, "site")
            label = getattr(scatterer, "label")
            u_iso = getattr(scatterer, "u_iso")
            if self._instructions:
                instruction = self._instructions.get(label, None)
            if instruction:
                x, y, z = map(float, numpy.dot(
                                  core.symop2numpy(core.xyzt2augmat(instruction).inv()), 
                                  numpy.array([[x0],[y0],[z0],[1.]])
                                  )[:3, 0])
            else:
                x, y, z = x0, y0, z0
            scatterer_intermediate = \
                    self._model_intermediate.get_by_label(label)
            if not scatterer_intermediate:
                print("Scatterer with the label {} is not found!".format(label))
                continue
            dx,dy,dz = [divmod(coor, self._modulo_n)[0] for coor in\
                                    getattr(scatterer_intermediate , "site")]
            if not scatterer_intermediate:
                raise KeyError("No scatterer with label: {} found !".format(label))
            #Index of the scatterer to change
            index = self._model_updated.get_index_by_label(label)
            setattr(self._model_updated.scatterers[index], "site", (x+dx,y+dy,z+dz)) 
            setattr(self._model_updated.scatterers[index], "u_iso", u_iso)
            diffarchitect.trans_xyz(
                                self._model_updated.scatterers[index],
                                core.symop2numpy(self._symop)
                                )

    def _save_model_updated(self):
        fbase, ext = self._model_path.split(".")
        model_updated_path = fbase + "_UPDATED.cif"
        print("Saving the updated model into: {}".format(model_updated_path))
        self._model_updated.save_CIF(model_updated_path)
                                  

    def update(self, sender):
        super(ModelUpdater, self).update(sender)
        self._update_model_updated()
        self._save_model_updated()
    
if __name__ == "__main__":
    #crystal_symmetry = crystal.symmetry(
    #            unit_cell = "4.2107 4.2107 13.093 90.0 90.0 90.0",
    #            space_group_symbol = "P 1")
    #fpath = "~/Dropbox/Python/Diffuse/disorder_1D/BEA.cif"
    #fr = ShelxModel(fpath)
    #print fr.unit_cell

    dirpath = "/Users/dima/Dropbox/Python/Diffuse/disorder_1D/BEA/fragment/refined_fragment"
    sender = SenderAdapter.tests(dirpath)    


    
