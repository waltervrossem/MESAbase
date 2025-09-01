#!/usr/bin/env python

import os
import numpy as np

from wsssss.inlists import create_grid as cg
from wsssss.inlists import inlists as inl

grid_name = 'grid'
grid = cg.MesaGrid(add_base_workdir=True, inlists_index=4)

mixture_inlist = 'inlist_a09'

grid.add_dir('src')
grid.add_file('history_columns.list')
grid.add_file('profile_columns.list')
grid.add_file('inlist_common')
grid.add_file(mixture_inlist)

grid.inlist['star_job']['read_extra_star_job_inlist(1)'] = True
grid.inlist['star_job']['extra_star_job_inlist_name(1)'] = 'inlist_common'
grid.inlist['star_job']['read_extra_star_job_inlist(2)'] = True
grid.inlist['star_job']['extra_star_job_inlist_name(2)'] = mixture_inlist

grid.inlist['controls']['read_extra_controls_inlist(1)'] = True
grid.inlist['controls']['extra_controls_inlist_name(1)'] = 'inlist_common'
grid.inlist['controls']['read_extra_controls_inlist(2)'] = True
grid.inlist['controls']['extra_controls_inlist_name(2)'] = mixture_inlist

grid.inlist['kap']['read_extra_kap_inlist(1)'] = True
grid.inlist['kap']['extra_kap_inlist_name(1)'] = 'inlist_common'
grid.inlist['kap']['read_extra_kap_inlist(2)'] = True
grid.inlist['kap']['extra_kap_inlist_name(2)'] = mixture_inlist

# Add Imbriani 14N(p,g)15O
grid.add_dir('rates_tables')
grid.add_file('rate_list.txt')

grid.kap['user_num_kap_Xs'] = 10
grid.kap['user_kap_Xs'] = 0.0e0, 0.1e0, 0.2e0, 0.35e0, 0.5e0, 0.7e0, 0.8e0, 0.9e0, 0.95e0, 1.0e0
grid.kap['user_num_kap_Zs'] = 13
grid.kap['user_kap_Zs'] = 0.000e0, 0.0001e0, 0.0003e0, 0.001e0, 0.002e0, 0.004e0, 0.01e0, 0.02e0, 0.03e0, 0.04e0, 0.06e0, 0.08e0, 0.100e0
grid.kap['user_num_kap_Xs_for_this_Z'] = 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 8

grid.kap['user_num_kap_CO_Xs'] = 5
grid.kap['user_kap_CO_Xs'] = 0.00e0, 0.03e0, 0.10e0, 0.35e0, 0.70e0
grid.kap['user_num_kap_CO_Zs'] = 8
grid.kap['user_kap_CO_Zs'] = 0.000e0, 0.001e0, 0.004e0, 0.010e0, 0.020e0, 0.030e0, 0.050e0, 0.100e0
grid.kap['user_num_kap_CO_Xs_for_this_Z'] = 5, 5, 5, 5, 5, 5, 5, 5

grid.kap['user_num_kap_lowT_Xs'] = 10
grid.kap['user_kap_lowT_Xs'] = 0.0e0, 0.1e0, 0.2e0, 0.35e0, 0.5e0, 0.7e0, 0.8e0, 0.9e0, 0.95e0, 1.0e0
grid.kap['user_num_kap_lowT_Zs'] = 13
grid.kap['user_kap_lowT_Zs'] = 0.000e0, 0.0001e0, 0.0003e0, 0.001e0, 0.002e0, 0.004e0, 0.01e0, 0.02e0, 0.03e0, 0.04e0, 0.06e0, 0.08e0, 0.100e0
grid.kap['user_num_kap_lowT_Xs_for_this_Z'] = 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 8
aFe_avail = np.linspace(-0.2, 0.6, 5)
aFe_avail_lowT = np.linspace(-0.2, 0.6, 9)


def inlist_finalize_function(unpacked_inlist):
    eos_file_prefix = 'mesa'
    aFe = unpacked_inlist['controls'][f'{non_mesa_key_start}_aFe']
    use_aFe = aFe_avail[np.argmin(np.abs(aFe_avail - aFe))]
    use_aFe_lowT = aFe_avail_lowT[np.argmin(np.abs(aFe_avail_lowT - aFe))]
    unpacked_inlist['kap']['kappa_file_prefix'] = f'agss09_afe_{use_aFe:+.1f}'
    unpacked_inlist['kap']['kappa_lowT_prefix'] = f'lowT_fa05_a09_aFe_{use_aFe_lowT:+.1f}'.replace('+0.', 'p').replace('-0.', 'm')
    unpacked_inlist['kap']['kappa_CO_prefix'] = f'agss09_afe_{use_aFe:+.1f}_co'
    return unpacked_inlist

# Add custom options here:
# e.g.
# grid.controls['initial_mass] = [1, 2]
grid.controls[f'{non_mesa_key_start}_aFe'] = 0.0  # [alpha/Fe] to use for grid

if __name__ == "__main__":
    grid.create_grid(grid_name)
    os.system(f'cp {__file__} {grid_name}/')
