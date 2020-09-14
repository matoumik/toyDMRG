#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:09:09 2020

@author: mikulas
"""

import DMRG as d

Dinfo = d.DMRG_info()
DMPS = d.MPS(Dinfo) 
DMPS.init_tensors()

DMPS.init_guess_rand(1.0)

print(DMPS.tensors)

