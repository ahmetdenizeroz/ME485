# -*- coding: utf-8 -*-
from utils.nb import dot
from utils.np import eps

import numpy as np
import re


def get_fsolver(name, be, cplargs):
    """
    docstring
    """
    fname = re.sub('\+', 'p', name)
    fname = re.sub('-', '_', fname)
    flux = eval('make_' + fname)(cplargs)
    return be.compile(flux)


def make_rusanov(cplargs):
    nvars     = cplargs['nfvars']
    gamma     = cplargs['gamma']
    flux_func = cplargs['flux_func']
    array     = cplargs['array']
    ndims     = cplargs['ndims']
    def fsolver(ul, ur, vl, vr, nf, fn):
        fl = array(nvars)
        fr = array(nvars)

        # this is u*phi*n 
        flux_func(ul, vl, nf, fl)
        flux_func(ur, vr, nf, fr)
        
        #---------------------------------#
        # TO DO: ????????????
        vnl = dot(vl, nf, ndims)
        vnr = dot(vr, nf, ndims)
        a = max(abs(vnl), abs(vnr))

        for i in range(nvars):
            fn[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * a * (ur[i] - ul[i])
        #---------------------------------#  

    return fsolver

def make_upwind(cplargs):
    nvars     = cplargs['nfvars']
    gamma     = cplargs['gamma']
    flux_func = cplargs['flux_func']
    array     = cplargs['array']
    ndims     = cplargs['ndims']
    def fsolver(ul, ur, vl, vr, nf, fn):
        print("upwindhear1")
        fl = array(nvars)
        print("upwindhear2")
        fr = array(nvars)
        print("upwindhear3")

        flux = dot(vl, nf, ndims)
        print("upwindhear6")

        flux_func(ul, vl, nf, fl)
        print("upwindhear4")

        flux_func(ur, vr, nf, fr)
        print("upwindhear5")

        # this is u*phi*n
        if flux > 0:
            fn[:] = fl[:]
            print("upwindhear7")
        else:
            fn[:] = fr[:]
            print("upwindhear8")
        
        #---------------------------------#  
        # complete the function
        #---------------------------------#  

    return fsolver