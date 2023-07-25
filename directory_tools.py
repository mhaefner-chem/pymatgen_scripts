# this function generates a name based on a structure
def directory_name_constructor(struc):
    from pymatgen.core import Composition
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    composition = Composition(struc.composition)

    # get the reduced sum formula
    formula = composition.reduced_formula
    formula = formula.replace('(','_')
    formula = formula.replace(')','_')

    # analyze space group
    space_group = SpacegroupAnalyzer(struc)
    HM_number = space_group.get_space_group_number()
    Z = composition.get_reduced_composition_and_factor()[1]

    # build name from formula, HM-number, and formula units in the unit cell
    return formula+"_SG"+str(HM_number)+"_Z"+str(Z)

# this function creates a directory
def make_directory(dir):
    import os
    if not os.path.exists(dir):
        os.makedirs(dir)
    else:
        print('The directory '+os.path.basename(dir)+' already exists and was skipped.')
    return

def find_head_directory():
    import os,sys
    # locate head directory
    while True:
        print(os.getcwd())
        if not os.path.isfile(".anchor"):
            os.chdir("..")
        else:
            print("Head directory found.")
            workdir = os.getcwd()
            break
        if os.getcwd() == "/":
            print("No anchor found.")
            sys.exit()
    return workdir