#!/usr/bin/env python

import os
import numpy as np

from wsssss.inlists import create_grid as cg
from wsssss.inlists import inlists as inl

grid = cg.MesaGrid(add_base_workdir=True, inlists_index=4)
grid.add_dir('src')
grid.add_file('history_columns.list')
grid.add_file('profile_columns.list')

mixture_inlist = 'inlist_a09'

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

# Add custom options here:


grid.create_grid('grid')
