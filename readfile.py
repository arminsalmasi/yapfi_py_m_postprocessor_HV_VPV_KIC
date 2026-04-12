import numpy as np
import os

def read_mapped_file(path, filename, dtype, optional=False):
    filepath = os.path.join(path, filename)
    if optional and not os.path.exists(filepath):
        return []
    with open(filepath, 'r') as f:
        if dtype is str:
            return f.read().split()
        return list(map(dtype, f))

def read_files(path):
    #path = os.getcwd() + '\\test2D'
    print(path)
# Read Files
    fin_volcentr_coord = read_mapped_file(path, 'FINITE_VOLUME_CENTROID_COORDINATES.TXT', float)
    chem_potentials = read_mapped_file(path, 'CHEMICAL_POTENTIALS.TXT', float)
    domain_size = read_mapped_file(path, 'DOMAIN_SIZE.TXT', float)
    grad_energy_contr = read_mapped_file(path, 'GRADIENT_ENERGY_CONTRIBUTION.TXT', float)
    mole_fractions = read_mapped_file(path, 'MOLE_FRACTIONS.TXT', float)
    n_elements = read_mapped_file(path, 'NUMBER_OF_ELEMENTS.TXT', int)
    n_gridpoints = read_mapped_file(path, 'NUMBER_OF_GRID_POINTS.TXT', int)
    n_phases = read_mapped_file(path, 'NUMBER_OF_PHASES.TXT', int)
    permeabilities = read_mapped_file(path, 'PERMEABILITIES.TXT', float)
    ph_field = read_mapped_file(path, 'PHASE_FIELD.TXT', float)
    ph_fractions = read_mapped_file(path, 'PHASE_FRACTIONS.TXT', float)
    time = read_mapped_file(path, 'TIME.TXT', float)
    el_names = read_mapped_file(path, 'ELEMENT_NAMES.TXT', str)
    ph_names = read_mapped_file(path, 'PHASE_NAMES.TXT', str)
    n_dimentions = read_mapped_file(path, 'DIMENSIONALITY.TXT', int)

    HCC = read_mapped_file(path, 'HCC_GPa.TXT', float, optional=True)
    K1C = read_mapped_file(path, 'K1C_MPa.TXT', float, optional=True)

 #Return Values
    return fin_volcentr_coord,  chem_potentials,  domain_size,  grad_energy_contr,  \
        mole_fractions,  n_elements,  n_gridpoints,  n_phases,  permeabilities,  \
        ph_field,  ph_fractions,  time,  el_names,  ph_names, n_dimentions, HCC, K1C