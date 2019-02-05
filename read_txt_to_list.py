import numpy as np
import os
from readfile import read_files
from loaddata import get_itemwise_scalars
from writevtk import write_vtk
from addbcs import add_scalars_limits
from addbcs import add_coords_limits
import re
##from tc_python import *

def main():
    #path = os.getcwd() + '\\8-2D-F285-TCFe-AIMD-FittedWithYapfi'
    #path = os.getcwd() + '\\LM-Yapfi-firstTry-20181221-files'
    path = os.getcwd() + '\\10-2D-F2275-TCFE-AIMD-FittedWithYapfi'
    #path = os.getcwd() + '\\8-1D-LM_2dGradSint-YAPFI-TCFE8-AIMD-1mm-YAPFIgrid-Compare_to_F285'

    [fin_volcentr_coord, chem_potentials, domain_size, grad_energy_contr, mole_fractions, n_elements,
     n_gridpoints, n_phases, permeabilities, ph_field, ph_fractions, time, el_names, ph_names, n_dimensions, hcc, k1c] = read_files(path)
    ntp = int(len(time))
    ngd = np.array(n_gridpoints)
    ngd_tot = np.prod(np.array(ngd))
    nel = n_elements[0]
    nph = n_phases[0]
    mole_fractions = np.array(mole_fractions)
    mole_fractions = np.reshape(mole_fractions, (ntp, ngd_tot*nel))
    array = np.array(chem_potentials)
    chem_potentials = array
    chem_potentials = np.reshape(chem_potentials, (ntp, ngd_tot*nel))
    ph_fractions = np.array(ph_fractions)
    ph_fractions = np.reshape(ph_fractions,(ntp, ngd_tot*nph))
    fin_volcentr_coord = np.array(fin_volcentr_coord)
    if len(hcc) != 0:
        hcc = np.reshape(hcc, (1,ngd_tot))
    if len(k1c) != 0:
        k1c = np.reshape(k1c, (1,ngd_tot))
    n_dim = n_dimensions[0]
    coords_limits, ngd_ternary = np.zeros(3), np.zeros(3, dtype=int)
    for i in range(n_dim):
        coords_limits[i] = domain_size[i]
        ngd_ternary[i] = ngd[i]
    if not os.path.exists(path + '\\BCvtk'):
        os.makedirs(path + '\\BCvtk')
    fin_volcentr_coord_with_limits = add_coords_limits(fin_volcentr_coord, ngd_ternary, coords_limits, n_dim)
    ngd_bc = np.array([], dtype=int)
    for i in range(n_dim):
        ngd_bc = np.append(ngd_bc, ngd[i]+2)
    for i in range(ntp):
        mf_tstp = get_itemwise_scalars(mole_fractions[i][:], nel, ngd_tot)
        chem_pot_tstp = get_itemwise_scalars(chem_potentials[i][:], nel, ngd_tot)
        phf_tstp = get_itemwise_scalars(ph_fractions[i][:], nph, ngd_tot)
        mf_bc = add_scalars_limits(mf_tstp[0, :], ngd)
        chem_pot_bc = add_scalars_limits(chem_pot_tstp[0, :], ngd)
        phf_bc = add_scalars_limits(phf_tstp[0, :], ngd)
        for elidx in range(1, nel):
            app_temp = []
            app_temp = add_scalars_limits(mf_tstp[elidx, :], ngd)
            mf_bc = np.vstack((mf_bc, app_temp))
            app_temp = []
            app_temp = add_scalars_limits(chem_pot_tstp[elidx, :], ngd)
            chem_pot_bc = np.vstack((chem_pot_bc, app_temp))
        for phidx in range(1, nph):
            app_temp = []
            app_temp = add_scalars_limits(phf_tstp[phidx, :], ngd)
            phf_bc = np.vstack((phf_bc, app_temp))
        header_str = write_vtk(fin_volcentr_coord_with_limits, ngd_bc, i, el_names, ph_names, mf_bc, chem_pot_bc, phf_bc)

        if (i == ntp-1) and (len(hcc) != 0):
            hcc_bc = add_scalars_limits(hcc[0][:], ngd)
            header_str = header_str + 'SCALARS ' + 'HCC' + ' Double 1' + '\n' + 'LOOKUP_TABLE default' + '\n' \
                   + re.sub('[\[\]]', '', np.array_str(hcc_bc[0][:])) + '\n'
        if (i == ntp-1) and (len(k1c) != 0):
            k1c_bc = add_scalars_limits(k1c[0][:], ngd)
            header_str = header_str + 'SCALARS ' + 'K1C' + ' Double 1' + '\n' + 'LOOKUP_TABLE default' + '\n' \
                         + re.sub('[\[\]]', '', np.array_str(k1c_bc[0][:])) + '\n'

        out_file_name = path + '\\BCvtk\\tstp_bc_' + str(i) + '.vtk'
        fout = open(out_file_name, "w")
        fout.write(header_str)
        fout.close()
0000
if __name__ == "__main__":
    main()







## without boundary limits
# header_str = write_vtk(fin_volcentr_coord, ngd, i, el_names, ph_names, mf_tstp, chem_pot_tstp, phf_tstp)
# if i == ntp:
#    HCC = np.array(HCC)
#    header_str = header_str + 'SCALARS ' + 'HCC' + ' Double 1' + '\n' + 'LOOKUP_TABLE default' + '\n' + re.sub('[\[\]]', '', np.array_str(HCC[0][:])) + '\n'
# out_file_name = path + '\\vtk\\tstp_' + str(i) + '.vtk'
# fout = open(out_file_name, "w")
# fout.write(header_str)
# fout.close()
# if not os.path.exists(path + '\\vtk'):
#   os.makedirs(path + '\\vtk')