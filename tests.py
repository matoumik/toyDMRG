#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:09:09 2020

@author: mikulas
"""

import DMRG as d

Dinfo = d.DMRG_info()

DMRG = d.DMRG(Dinfo)
print(DMRG.MPS.M)
DMRG.MPS.init_guess_rand(1.0)
#DMRG.solve_onesite(0)
#DMRG.solve_onesite(1)
#DMRG.solve_onesite(0)
#DMRG.solve_onesite(1)
#DMRG.solve_onesite(0)
#DMRG.solve_onesite(15)
#DMRG.solve_onesite(14)
#DMRG.sweep_backw()
DMRG.do_sweep()
for i in range(10):
    DMRG.do_sweep()


#print(DMPS.tensors)

