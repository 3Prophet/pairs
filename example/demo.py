###############################################################################
# This module demonstrates basic functionality of the API                     #  
###############################################################################

import os

import pairs.corrtable as corrtable
from pairs.symmlib import core
from pairs.stacking import diffarchitect 

def stackingDemo():
    #--------------------------------------------------------------
    #Average layer group
    L = map(core.modulo_n,
            core.applySymops(['x,y,z', 'x+1/3,y,z','x+2/3,y,z'])(
            core.applySymops(['x,y,z', 'x,y+1/3,z','x,y+2/3,z'])(
            core.shiftOrigin('x,y,z+1/8')(core.symopSymList(115)))))

    ##--------------------------------------------------------------
    fragment = corrtable.Fragment(stab = core.shiftOrigin('x,y,z+1/8')(
                                        #['x,y,z', '-y,-x,-z']
                                        core.symopSymList(115)
                                        ))

    layer_0 = corrtable.Layer(stab = L, fragments = [fragment])
    layer_1 = corrtable.Layer.leftmult(layer_0, "x,y,-z+1/2")

    bilayer_01 = corrtable.Bilayer([layer_0, layer_1])

    table = corrtable.PairCorrTab(bilayer_01)
    table.corr_analysis()
    print table
    table.fill_with_values({
                    1: '(1-2*k-2*l-4*m)*1/9',
                    2: 'k*1/9',
                    3: 'l*1/9',
                    4: 'm*1/9'
                    })
    print table
    return bilayer_01 

def architectDemo(bilayer):
    dirPath = os.getcwd()
    fPath = os.path.join(dirPath, "BEA.cif")
    demoDirPath = os.path.join(dirPath, "examplefiles") #dir to store demos
    fragPathOut = os.path.join(demoDirPath, "BEA_fragment.cif")
    layerPathOut = os.path.join(demoDirPath, "BEA_layer_0.cif")
    bilayerPathOut = os.path.join(demoDirPath, "BEA_bilayer_01.cif")
    pairPathOut = os.path.join(demoDirPath, "BEA_pair_12.cif")
    #--------------------------------------------------------------
    #Cutting the fragment
    initializer = diffarchitect.Initializer(fPath, expand_to_p1 = True)
    fragment = initializer.generate_fragment().limit_to_range(
            x = (0, .999), y = (0, .999), z = (0, 0.25))
    #Generating unique labels
    fragment.gen_unique_labels()
    #Shifting the fragment
    fragment.transform('x-1/3,y-2/3,z')
    fragment.transform_cond('x-1,y,z', 'x>0.51')
    fragment.transform_cond('x,y-1,z', 'y>0.49')

    model_f = diffarchitect.Model(initializer.symmetry, [fragment])
    model_f.save_CIF(fragPathOut)
    #return model_f

    #--------------------------------------------------------------
    #Generating layer:
    fragments_0 = map(fragment.clone,
            [frag.coset[0] for frag in bilayer.layers[0].fragments])
    layer_0 = diffarchitect.Model(
        symmetry = initializer.symmetry,
        fragments = fragments_0
        )

    layer_0.save_CIF(layerPathOut)
    #--------------------------------------------------------------
    #Generating bilayer:
    fragments_1 = map(fragment.clone,
            [frag.coset[0] for frag in bilayer.layers[1].fragments])
    bilayer_01 = diffarchitect.Model(
            symmetry = initializer.symmetry,
            fragments = core.flatMap([fragments_0, fragments_1])
            )
    bilayer_01.save_CIF(bilayerPathOut)
    #---------------------------------------------------------------
    #Generating a pair
    pair_fragments = map(fragment.clone, [
                                    bilayer.layers[0].fragments[1].coset[0],
                                    bilayer.layers[1].fragments[2].coset[0]
                                    ])
    pair_12 = diffarchitect.Model(
            symmetry = initializer.symmetry,
            fragments = pair_fragments
            )
    pair_12.save_CIF(pairPathOut)

    pair_23 = diffarchitect.Model(
            symmetry = initializer.symmetry,
            fragments = map(fragment.clone, ['x+1/3,y,z',  'x,y,-z+1/2'])
            )
    pair_23.save_CIF(pairPathOut)    

if __name__ == "__main__":
    stackingDemo()
    architectDemo(stackingDemo())



    
