import numpy as np
import os

def read_files(path):
    #path = os.getcwd() + '\\test2D'
    print(path)
# Read Files
    fname = ''.join(path + '/FINITE_VOLUME_CENTROID_COORDINATES.TXT')
    with open(fname, 'r') as f:
        fin_volcentr_coord = map(float, f)
    f.close()
    fname = ''.join(path + '/CHEMICAL_POTENTIALS.TXT')
    with open(fname, 'r') as f:
        chem_potentials = map(float, f)
    f.close()
    fname = ''.join(path + '/DOMAIN_SIZE.TXT')
    with open(fname, 'r') as f:
        domain_size = map(float, f)
    f.close()
    fname = ''.join(path + '/GRADIENT_ENERGY_CONTRIBUTION.TXT')
    with open(fname, 'r') as f:
        grad_energy_contr = map(float, f)
    f.close()
    fname = ''.join(path + '/MOLE_FRACTIONS.TXT')
    with open(fname, 'r') as f:
        mole_fractions = map(float, f)
    f.close()
    fname = ''.join(path + '/NUMBER_OF_ELEMENTS.TXT')
    with open(fname, 'r') as f:
        n_elements = map(int, f)
    f.close()
    fname = ''.join(path + '/NUMBER_OF_GRID_POINTS.TXT')
    with open(fname, 'r') as f:
        n_gridpoints = map(int, f)
    f.close()
    fname = ''.join(path + '/NUMBER_OF_PHASES.TXT')
    with open(fname, 'r') as f:
        n_phases = map(int, f)
    f.close()
    fname = ''.join(path + '/PERMEABILITIES.TXT')
    with open(fname, 'r') as f:
        permeabilities = map(float, f)
    f.close()
    fname = ''.join(path + '/PHASE_FIELD.TXT')
    with open(fname, 'r') as f:
        ph_field = map(float, f)
    f.close()
    fname = ''.join(path + '/PHASE_FRACTIONS.TXT')
    with open(fname, 'r') as f:
        ph_fractions = map(float, f)
    f.close()
    fname = ''.join(path + '/TIME.TXT')
    with open(fname, 'r') as f:
        time = map(float, f)
    f.close()	
    fname = ''.join(path + '/ELEMENT_NAMES.TXT')
    with open(fname, 'r') as f:
        el_names = f.read().split()
    f.close()	
    fname = ''.join(path + '/PHASE_NAMES.TXT')
    with open(fname, 'r') as f:
        ph_names = f.read().split()
    f.close()
    fname = ''.join(path + '/DIMENSIONALITY.TXT')
    with open(fname, 'r') as f:
        n_dimentions = map(int, f)
    f.close()

    HCC = []
    K1C = []
    HCCExist = os.path.exists(path + '/HCC_GPa.TXT')
    K1CExist = os.path.exists(path + '/K1C_MPa.TXT')
    if HCCExist:
        fname = ''.join(path + '/HCC_GPa.TXT')
        with open(fname, 'r') as f:
            HCC = map(float, f)
        f.close()
    if K1CExist:
        fname = ''.join(path + '/K1C_MPa.TXT')
        with open(fname, 'r') as f:
            K1C = map(float, f)
        f.close()

 #Return Values
    return fin_volcentr_coord,  chem_potentials,  domain_size,  grad_energy_contr,  \
        mole_fractions,  n_elements,  n_gridpoints,  n_phases,  permeabilities,  \
        ph_field,  ph_fractions,  time,  el_names,  ph_names, n_dimentions, HCC, K1C