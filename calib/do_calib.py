#!/usr/bin/env python

import os
import sys
import time

import numpy as np
from scipy import optimize as op
from scipy import interpolate as ip

from wsssss.inlists import create_grid as cg
from wsssss.inlists import inlists as inl
from wsssss import load_data as ld
from wsssss._bin.check_grid import get_mesa_termcode

if not 'SLURM_JOB_ID' in os.environ.keys():
    os.environ['SLURM_JOB_ID'] = time.strftime("%Y-%m-%d_%H%M%S")

args = sys.argv[1:]
solar_mixture = args[0]
if len(args) == 2:
    version = args[1]
else:
    version = '24081'

defaults = os.path.abspath(f'../{version}')

if solar_mixture not in ['agss09', 'mb22', 'gs98']:
    raise ValueError('Solar mixture not agss09 or mb22.')
print(f'{solar_mixture=}')
print(f'{version=}')

_kap_file = {'agss09': 'a09',
             'mb22': 'oplib_mb22',
             'gs98': 'gs98'}
_kap_lowT = {'agss09': 'a09p',
             'mb22': 'mb22',
             'gs98': 'gs98'}
_kap_CO = {'agss09': 'a09',
           # 'mb22': 'mb22'}
           'mb22': 'gs98',  # No mb22 CO tables, gs98 is close to mb22 z/x (0.0231 vs 0.0225)
           'gs98': 'gs98'}
_initial_zfracs = {'agss09': 6,
                   'mb22': 9,
                   'gs98': 3}

_ZXsol = {'agss09': 0.0178,  # http://adsabs.harvard.edu/abs/2009ApJ...705L.123S
          'mb22': 0.0225,
          'gs98': 0.0229
          }
_e_ZXsol = {'agss09':0.0011,  # For agss09 use same fractional error as mb22
            'gs98': 0.0023,  # 10% error from gs98 paper
            'mb22': 0.0014}

if not 'SLURM_JOB_ID' in os.environ.keys():
    os.environ['SLURM_JOB_ID'] = time.strftime("%Y-%m-%d_%H%M%S")

grid_name = f'solar_calib_{solar_mixture}_{version}_{os.environ["SLURM_JOB_ID"]}'
grid = cg.MesaGrid(add_base_workdir=True, inlists_index=4)
grid.add_dir(f'{defaults}/src')
grid.add_file(f'{defaults}/history_columns.list')
grid.add_file(f'{defaults}/profile_columns.list')

grid.inlist['star_job']['read_extra_star_job_inlist(1)'] = True
grid.inlist['star_job']['extra_star_job_inlist_name(1)'] = f'{defaults}/inlist_common'

grid.inlist['controls']['read_extra_controls_inlist(1)'] = True
grid.inlist['controls']['extra_controls_inlist_name(1)'] = f'{defaults}/inlist_common'

grid.inlist['kap']['read_extra_kap_inlist(1)'] = True
grid.inlist['kap']['extra_kap_inlist_name(1)'] = f'{defaults}/inlist_common'

grid.inlist['star_job']['read_extra_star_job_inlist(5)'] = True
grid.inlist['star_job']['extra_star_job_inlist_name(5)'] = 'inlist_calib'

grid.inlist['controls']['read_extra_controls_inlist(5)'] = True
grid.inlist['controls']['extra_controls_inlist_name(5)'] = 'inlist_calib'

grid.inlist['kap']['read_extra_kap_inlist(5)'] = True
grid.inlist['kap']['extra_kap_inlist_name(5)'] = 'inlist_calib'

# Add Imbriani 14N(p,g)15O
grid.add_dir(f'{defaults}/rates_tables')
grid.add_file(f'{defaults}/rate_list.txt')

# Add custom options here:
grid.star_job['initial_zfracs'] = _initial_zfracs[solar_mixture]
grid.kap['kap_file_prefix'] = f'{_kap_file[solar_mixture]}'
grid.kap['kap_lowT_prefix'] = f'lowT_fa05_{_kap_lowT[solar_mixture]}'
grid.kap['kap_CO_prefix'] = f'{_kap_CO[solar_mixture]}_co'

try:
    grid.create_grid(grid_name)
except FileNotFoundError:  # It can't find inlist_calib but we create it when running the calibration.
    pass

os.system(f'cp {__file__} {grid_name}/')

run_dir = os.path.abspath(f'{grid_name}/0000')
os.chdir(run_dir)

ts = time.time()

max_age = 5.5e9
max_age_inlist_str = f'{max_age:.1e}'.lower().replace('e', 'd').replace('+', '')
min_age = 4.e9
# target
age_target = 4567.3e6  # https://www.science.org/doi/10.1126/science.1226919
precision_age = 0.16e6

logL_target = 0
precision_logL = 0.0014/3.8275 * np.log10(np.e)  # IAU 2015 Resolution B3 references

logR_target = 0
precision_logR = 140/695658 * np.log10(np.e)  # IAU 2015 Resolution B3 references

Z_X = _ZXsol[solar_mixture]
precision_Z_X =_e_ZXsol[solar_mixture]

# Basu 2004
Y_surf = 0.2485
precision_Y = 0.0034

target = np.array([logL_target, logR_target, Z_X, Y_surf, age_target])
precision = np.array([precision_logL, precision_logR, precision_Z_X, precision_Y, precision_age])

# inputs
a_MLT_ini = 2.0
Y_ini = 0.273
Z_ini = 0.8 * Z_X

x0 = np.array([a_MLT_ini, Y_ini, Z_ini])
x_scale = x0
calib_fname = os.path.abspath(f'{run_dir}/../calib_{solar_mixture}.dat')
if os.path.exists(calib_fname):
    os.remove(calib_fname)
print(calib_fname)

def do_run(x, _counter=[0]):
    a_MLT_ini, Y_ini, Z_ini = x
    print(f'Starting =======')
    inlist_str = f'''
&star_job
   save_photo_when_terminate = .true.
/

&controls
   initial_z = {Z_ini}
   initial_y = {Y_ini}
   mixing_length_alpha = {a_MLT_ini}
   max_age = {max_age_inlist_str}
   num_adjusted_dt_steps_before_max_age = 0
   star_history_name = 'history_{_counter[0]:03}.data'
/

&kap
   Zbase = {Z_ini}
/
'''
    with open(f'inlist_calib', 'w') as handle:
        handle.write(inlist_str)
        os.fsync(handle)  # write to file before doing anything else
    ierr = os.system(f'./rn > tmp_out')
    if ierr != 0:
        raise RuntimeError(ierr)

    tmp = lambda : None
    tmp.grid_dir = '../'
    tmp.out_file = f'tmp_out'
    term_code = get_mesa_termcode('0000', tmp)[1]
    if term_code != 'max_age':  # Something went wrong, try again with different mesh_delta
        inlist_str.replace('/\n\n&kap', 'mesh_delta_coeff = 0.98d0/\n\n&kap')
        with open(f'inlist_calib', 'w') as handle:
            handle.write(inlist_str)
            os.fsync(handle)  # write to file before doing anything else
        ierr = os.system(f'./rn > tmp_out')
        term_code = get_mesa_termcode('0000', tmp)[1]
        if term_code != 'max_age':
            raise ValueError(f'Incorrect termination code: {term_code}')

    hist = ld.History(f'LOGS/history_{_counter[0]:03}.data', verbose=False)
    os.system(f'cp photos/x{hist.data.model_number[-1]%1000:05} photos/final{_counter[0]:03}')
    mask = hist.data.star_age >= min_age
    mask[max(0, np.where(mask)[0][0] - 1)] = True  # Add one point before min_age to avoid bounds_error in interpolation

    X_surf = hist.data.surface_h1[mask] + hist.data.surface_h2[mask]
    Y_surf = hist.data.surface_he3[mask] + hist.data.surface_he4[mask]
    Z_surf = 1 - X_surf - Y_surf
    ZX_surf = Z_surf / X_surf

    logL = np.log10(hist.data.luminosity[mask])
    logR = np.log10(hist.data.photosphere_r[mask])
    Teff = hist.data.effective_T[mask]

    age = hist.data.star_age[mask]

    weighted_resi = (np.array([logL, logR, ZX_surf, Y_surf, age]) - target[:, np.newaxis]) / precision[:, np.newaxis]
    ip_res = ip.interp1d(age, weighted_resi, kind='quadratic')
    ip_fun = lambda x: ip_res(x[0])  # Get correct ouptut shape for least_squares
    res = op.least_squares(ip_fun, x0=4.65e9, bounds=[min_age, max_age])  # fit for best age
    best_age = res['x'][0]

    ip_vals = ip.interp1d(age, np.array([logL, logR, X_surf, Y_surf, ZX_surf, Teff]), kind='quadratic')
    logL, logR, X_surf, Y_surf, ZX_surf, Teff = ip_vals(best_age)
    cost = res.cost

    print(f'Cost     = {cost: e}')
    print(f'a_mlt    = {a_MLT_ini: f}')
    print(f'X_init   = {1 - Y_ini - Z_ini: f}')
    print(f'Y_init   = {Y_ini: f}')
    print(f'Z_init   = {Z_ini: f}')
    print(f'Best age = {best_age: e}')
    print(f'logL     = {logL: e}')
    print(f'logR     = {logR: e}')
    print(f'T_eff    = {Teff: f}')
    print(f'X        = {X_surf: f}')
    print(f'Y        = {Y_surf: f}')
    print(f'Z/X      = {ZX_surf: f}')

    if os.path.exists(calib_fname):
        file_mode = 'a'
    else:
        file_mode = 'w'
    column_names = 'Run, Cost, logL, logR, a_mlt, Y_ini, Z_ini, Teff, X_surf, Y_surf, ZX_surf, Age\n'
    with open(calib_fname, file_mode) as handle:
        if file_mode == 'w':  # Also write header
            handle.write(column_names)
        handle.write(','.join([f'{int(_counter[0]): 3}, {cost: e}', f'{logL: e}', f'{logR: e}',
                               f'{a_MLT_ini: f}', f'{Y_ini: f}', f'{Z_ini: f}',
                               f'{Teff: f}', f'{X_surf: f}', f'{Y_surf: f}', f'{ZX_surf: f}', f'{best_age: e}']) + '\n')
    _counter[0] += 1
    return res.fun

os.system('./clean && ./mk')

res = op.least_squares(do_run, x0, ftol=1e-3, x_scale=x_scale, diff_step=0.005)

os.system(f'cp {calib_fname} {calib_fname}.bak')

calib_dat = np.genfromtxt(calib_fname, delimiter=',', names=True)
columns = list(calib_dat.dtype.fields.keys())
order = np.argsort(calib_dat['Cost'])

with open(calib_fname, 'w') as handle:
    handle.write(', '.join(columns) + '\n')
    for i, line in enumerate(calib_dat[order]):
        run, cost, logL, logR, a_mlt, Y_ini, Z_ini, Teff, X_surf, Y_surf, ZX_surf, best_age = line
        run = int(run)
        X_ini = 1 - Z_ini - Y_ini
        run_idx = order[i]
        handle.write(','.join([f'{run: 3}', f'{cost: e}', f'{logL: e}', f'{logR: e}',
                               f'{a_mlt: f}', f'{Y_ini: f}', f'{Z_ini: f}',
                               f'{Teff: f}', f'{X_surf: f}', f'{Y_surf: f}', f'{ZX_surf: f}', f'{best_age: e}']) + '\n')

        if i == 0:
            print('========================')
            print('Best Fit:')
            print(f'Run      = {run: 3}')
            print(f'Cost     = {cost: e}')
            print(f'Best age = {best_age: e}')
            print(f'a_mlt    = {a_mlt: f}')
            print(f'X_init   = {1 - Y_ini - Z_ini: f}')
            print(f'Y_init   = {Y_ini: f}')
            print(f'Z_init   = {Z_ini: f}')
            print(f'Z/X_init = {Z_ini/X_ini: f}')
            print(f'logL     = {logL: e}')
            print(f'logR     = {logR: e}')
            print(f'T_eff    = {Teff: f}')
            print(f'X        = {X_surf: f}')
            print(f'Y        = {Y_surf: f}')
            print(f'Z        = {X_surf * ZX_surf: f}')
            print(f'Z/X      = {ZX_surf: f}')

# Write mixture specific inlist
line = calib_dat[order][0]
run, cost, logL, logR, a_mlt, Y_ini, Z_ini, Teff, X_surf, Y_surf, ZX_surf, best_age = line
run = int(run)
X_ini = 1 - Z_ini - Y_ini
run_idx = order[i]

s = (f'! {grid_name}\n'
f'! {column_names}\n'
f'! ')
s += ','.join([f'{run: 3}', f'{cost: e}', f'{logL: e}', f'{logR: e}', f'{a_mlt: f}', f'{Y_ini: f}', f'{Z_ini: f}',
                       f'{Teff: f}', f'{X_surf: f}', f'{Y_surf: f}', f'{ZX_surf: f}', f'{best_age: e}']) + '\n'
s += '\n'
s+= f'''
&star_job
    initial_zfracs = {_initial_zfracs[solar_mixture]} ! {solar_mixture}
    change_net = .true.
    new_net_name = 'pp_and_cno_extras.net'
/ ! end of star_job namelist

&kap
    kap_file_prefix = '{_kap_file[solar_mixture]}'
    kap_lowT_prefix = 'lowT_fa05_{_kap_lowT[solar_mixture]}'
    kap_CO_prefix = '{_kap_CO[solar_mixture]}_co'
    Zbase = {Z_ini}

/ ! end of kap namelist

&controls
    mixing_length_alpha = {a_mlt}
    initial_y = {Y_ini}
    initial_z = {Z_ini}
/ ! end of controls namelist
'''
inlist_name = 'inlist_' + {'agss09': 'a09',
               'gs98': 'gs98',
               'mb22': 'mb22'}[solar_mixture]
with open(os.path.abspath(f'{run_dir}/../{inlist_name}'), 'w') as handle:
    handle.write(s)

# Do full run at end
Z_ini = calib_dat[order][0][6]
Y_ini = calib_dat[order][0][5]
a_MLT_ini = calib_dat[order][0][4]
inlist_str = f'''
&star_job
   write_profile_when_terminate = .true.
   filename_for_profile_when_terminate = 'final.data'
   save_pulse_data_when_terminate = .true.
   save_pulse_data_filename = 'final.data.GYRE'
/

&controls
   initial_z = {Z_ini}
   initial_y = {Y_ini}
   mixing_length_alpha = {a_MLT_ini}
   max_age = 1d36
   photosphere_r_upper_limit = 20
   num_adjusted_dt_steps_before_max_age = 0
   star_history_name = 'history.data'
   pulse_data_format = 'GYRE'
   write_profiles_flag = .true.
/

&kap
   Zbase = {Z_ini}
/
'''
with open(f'inlist_calib', 'w') as handle:
    handle.write(inlist_str)
    os.fsync(handle)  # write to file before doing anything else
ierr = os.system(f'./rn > out_{solar_mixture}')
tt = time.time() - ts
print(f'TotalTime= {int(tt//3600)}h{(int(tt%3600))//60}m{int(tt%60)}s')
