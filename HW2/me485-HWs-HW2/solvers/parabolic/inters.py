# -*- coding: utf-8 -*-
from numba.np.arrayobj import np_shape
from numpy.f2py.auxfuncs import throw_error

from solvers.base import BaseIntInters, BaseBCInters, BaseMPIInters
from backends.types import Kernel, NullKernel
from solvers.parabolic.bcs import get_bc

from utils.np import npeval
from utils.nb import dot
import functools as fc
import numpy as np
import re

#-------------------------------------------------------------------------------#    
class ParabolicIntInters(BaseIntInters):
    def construct_kernels(self, elemap, impl_op):
        # View of elemenet array
        self._fpts = fpts = [cell.fpts for cell in elemap.values()]
        dfpts = [cell.grad for cell in elemap.values()]
        nele = len(fpts)

        self._correction = self.cfg.get('solver', 'correction', 'minimum')
       
        # Array for gradient at face
        self._gradf  = gradf   = np.empty((self.ndims, self.nvars, self.nfpts))
        self._weight = weight  = np.empty((self.nfpts))

        if self.order > 1:
            self.compute_delu = Kernel(self._make_delu(), *fpts)
            dxf = [cell.dxf for cell in elemap.values()]        
            self.compute_weight(dxf)
        else:
            self.compute_delu = NullKernel

        # Kernel to compute gradient at face (Averaging gradient)
        self.compute_grad_at_face = Kernel(
                                    self._make_grad_at_face(nele), gradf, *fpts, *dfpts)

        muf = np.empty(self.nfpts)
        self.compute_flux = Kernel(self._make_flux(nele), muf, gradf, *fpts)


#-------------------------------------------------------------------------------#    
    def compute_weight(self, dxf):
        nface, ndims = self.nfpts, self.ndims
        lt, le, lf = self._lidx
        rt, re, rf = self._ridx
        for idx in range(self.nfpts):
            lti, lfi, lei = lt[idx], lf[idx], le[idx]
            rti, rfi, rei = rt[idx], rf[idx], re[idx]

            dl = 0.0
            dr = 0.0
            for jdx in range(ndims):
                xjl = dxf[lti][lfi, jdx, lei]
                xjr = dxf[rti][rfi, jdx, rei]
                dl += xjl*xjl
                dr += xjr*xjr

            dl  = np.sqrt(dl)
            dr  = np.sqrt(dr)

            self._weight[idx] = dr/(dl+dr)

#-------------------------------------------------------------------------------#    
    def _make_flux(self, nele):
        ndims, nfvars = self.ndims, self.nfvars
        lt, le, lf = self._lidx
        rt, re, rf = self._ridx
        nf, sf = self._vec_snorm, self._mag_snorm

        # Inverse distance between the elements
        inv_ef    = self._rcp_dx
        
        # Unit vector connecting cell centers 
        ef = self._dx_adj * inv_ef

        # Correction method to be used
        correction = self._correction

        # Compiler arguments
        array = self.be.local_array()
        
        # Get compiled function of viscosity and viscous flux
        compute_mu = self.ele0.mu_container()

        def comm_flux(i_begin, i_end, muf, gradf, *uf):
            # Parse element views (fpts, grad)
            du    = uf[:nele]
            #Sf = nf * sf
            for idx in range(i_begin, i_end):
                muf = compute_mu()
                lti, lfi, lei = lt[idx], lf[idx], le[idx]
                rti, rfi, rei = rt[idx], rf[idx], re[idx]
                #fn = array(nfvars)
                fn = 0
                Sf = nf[:, idx] * sf[idx]
                # Ef sanırım burda hesaplanacak
                if correction == "minimum":
                    # Ef = np.dot(ee, Sfe) * ee
                    #Ef = dot(ef[:, idx], Sf[:, idx], ndims) #* ef[:, idx]
                    Ef = dot(ef[:, idx], Sf, ndims)  # * ef[:, idx]
                elif correction == "orthogonal":
                    Ef = sf[idx] #* ef[:, idx]
                elif correction == "over_relaxed":
                    # Ef = (sf[idx] * sf[idx]) / (np.dot(ee, nfe)) * ee
                    Ef = sf[idx] * (dot(nf[:, idx], nf[:, idx], ndims) / dot(ef[:, idx], nf[:, idx], ndims)) #* ef[:, idx]
                else:
                    print("Wrong correction type")
                #Tfe = Sf[:, idx] - Ef * ef[:,idx]
                Tfe = Sf - Ef * ef[:, idx]
                #print("tfe", Tfe)
                #temporary = np.linalg.norm(Ef) * ((du[lti][lfi][0][lei]-du[rti][lfi][0][rei])*inv_ef[idx]) + dot(gradf[:, 0, lei], Tfe, ndims)
                #temporary = np.linalg.norm(Ef) * (du[lti][lfi][0][lei] * inv_ef[idx]) + dot(gradf[:, 0, lei], Tfe, ndims)
                #temporary = Ef * (du[lti][lfi][0][lei] * inv_ef[idx]) + dot(gradf[:, 0, idx], Tfe, ndims)
                fn = -1 * muf * Ef * (du[lti][lfi][0][lei] * inv_ef[idx]) + dot(gradf[:, 0, idx], Tfe, ndims)
                #print("\n", gradf[:, 0, idx])
                #print("fn", fn)
                # Daha sonra burda flux hesaplanmalı
                uf[lti][lfi, 0, lei] =  fn
                uf[rti][rfi, 0, rei] = -fn
                #print("int flux final uf", uf)

        return self.be.make_loop(self.nfpts, comm_flux)


#-------------------------------------------------------------------------------#    
    def _make_grad_at_face(self, nele):
        nvars, ndims = self.nvars, self.ndims
        lt, le, lf = self._lidx
        rt, re, rf = self._ridx    
        # Inverse distance between the cell center
        weight    = self._weight
        # Stack-allocated array
        array = self.be.local_array()

        def grad_at(i_begin, i_end, gradf, *uf):
            # Parse element views (fpts, grad)
            du    = uf[:nele]
            gradu = uf[nele:]
            #print("uf", uf)
            #print("du", du)
            #print("gradu int", gradu)
            #print("lt", lt, np.shape(lt))
            #print("rt", rt, np.shape(rt))
            #print("le", le, np.shape(le))
            #print("re", re, np.shape(re))
            #print("lf", lf, np.shape(lf))
            #print("rf", rf, np.shape(rf))
            #print("weight", weight, np.shape(weight))
            #print("gradf int", gradf, np.shape(gradf))
            for idx in range(i_begin, i_end):
                #*************************#
                lti, lfi, lei = lt[idx], lf[idx], le[idx]
                rti, rfi, rei = rt[idx], rf[idx], re[idx]
                wl, wr = weight[idx], (1 - weight[idx])
                for dimension in range(ndims):
                    for variable in range(nvars):
                        ql = gradu[lti][dimension][variable][lei]
                        qr = gradu[rti][dimension][variable][rei]
                        gradf[dimension][variable][idx] = wl*ql + wr*qr
                #*************************#
            #print("gradf", gradf, np.shape(gradf))
        return self.be.make_loop(self.nfpts, grad_at)

#-------------------------------------------------------------------------------#    
    def _make_delu(self):
        nvars = self.nvars
        lt, le, lf = self._lidx
        rt, re, rf = self._ridx

        def compute_delu(i_begin, i_end, *uf):
            for idx in range(i_begin, i_end):
                lti, lfi, lei = lt[idx], lf[idx], le[idx]
                rti, rfi, rei = rt[idx], rf[idx], re[idx]

                for jdx in range(nvars):
                    ul = uf[lti][lfi, jdx, lei]
                    ur = uf[rti][rfi, jdx, rei]
                    du = ur - ul
                    uf[lti][lfi, jdx, lei] =  du
                    uf[rti][rfi, jdx, rei] = -du

        return self.be.make_loop(self.nfpts, compute_delu)

#-------------------------------------------------------------------------------#    
class ParabolicBCInters(BaseBCInters):
    _get_bc = get_bc
    def construct_bc(self):
        # Parse BC function name
        bcf = re.sub('-', '_', self.name)

        # Constants for BC function
        if self._reqs:
            bcsect = 'soln-bcs-{}'.format(self.bctype)
            bcc = {k: npeval(self.cfg.getexpr(bcsect, k, self._const))
                   for k in self._reqs}
        else:
            bcc = {}

        bcc['ndims'], bcc['nvars'], bcc['nfvars'] = self.ndims, self.nvars, self.nfvars



        bcc.update(self._const)
        
        # Get bc from `bcs.py` and compile them
        self.bc = self._get_bc(self.be, bcf, bcc)

#-------------------------------------------------------------------------------#    
    def construct_kernels(self, elemap, impl_op):
        self.construct_bc()

        self._correction = self.cfg.get('solver', 'correction', 'minimum')
        # View of elemenet array
        self._fpts = fpts = [cell.fpts for cell in elemap.values()]
        dfpts = [cell.grad for cell in elemap.values()]
        nele = len(fpts)

       # Gradient at face
        self._gradf = gradf = np.empty((self.ndims, self.nvars, self.nfpts))

        # Kernel to compute differnce of solution at face
        self.compute_delu = Kernel(self._make_delu(), *fpts)

        # Kernel to compute gradient at face (Averaging gradient)
        self.compute_grad_at_face = Kernel(
            self._make_grad_at_face(nele), gradf, *fpts, *dfpts
        )
  
        # Save viscosity on face (for implicit operator)
        muf = np.empty(self.nfpts)

        # Kernel to compute flux
        self.compute_flux = Kernel(self._make_flux(nele), muf, gradf, *fpts)

#-------------------------------------------------------------------------------#    
    def _make_flux(self, nele):
        ndims, nfvars = self.ndims, self.nfvars
        lt, le, lf = self._lidx
        nf, sf = self._vec_snorm, self._mag_snorm

        # Mangitude and direction of the connecting vector
        inv_ef = self._rcp_dx
        ef = self._dx_adj * inv_ef
        # avec = self._vec_snorm/np.einsum('ij,ij->j', ef, self._vec_snorm)

        correction = self._correction

        # Compiler arguments
        array = self.be.local_array()
       
        # Get compiled function of viscosity and viscous flux
        compute_mu = self.ele0.mu_container()
        # Get bc function 
        bc = self.bc

        def comm_flux(i_begin, i_end, muf, gradf, *uf):
            # Parse element views (fpts, grad)
            du    = uf[:nele]
            #Sf = nf * sf
            for idx in range(i_begin, i_end):
                muf = compute_mu()
                lti, lfi, lei = lt[idx], lf[idx], le[idx]
                #fn = np.empty(nfvars)
                fn = 0
                Sf = nf[:, idx] * sf[idx]
                # Ef sanırım burda hesaplanacak
                if correction == "minimum":
                    #Ef = np.dot(ee, Sfe) * ee
                    #Ef = dot(ef[:,idx], Sf[:,idx], ndims) #* ef[:,idx]
                    Ef = dot(ef[:, idx], Sf, ndims)  # * ef[:,idx]
                elif correction == "orthogonal":
                    Ef = sf[idx] #* ef[:,idx]
                elif correction == "over_relaxed":
                    #Ef = (sf[idx] * sf[idx]) / (np.dot(ee, nfe)) * ee
                    Ef = sf[idx] * (dot(nf[:,idx], nf[:,idx], ndims)/dot(ef[:,idx], nf[:,idx], ndims)) #* ef[:,idx]
                else:
                    print("Wrong correction type")
                #Tfe = Sf[:, idx] - Ef * ef[:,idx]
                Tfe = Sf - Ef * ef[:, idx]
                #temporary = np.linalg.norm(Ef) * (du[lti][lfi][0][lei] * inv_ef[idx]) + dot(gradf[:,0,lfi],Tfe, ndims)
                temporary = Ef * (du[lti][lfi][0][lei] * inv_ef[idx]) + dot(gradf[:, 0, idx], Tfe, ndims)
                fn = -1 * muf * temporary
                # Daha sonra burda flux hesaplanmalı
                uf[lti][lfi, 0, lei] = fn
                # print("int flux final uf", uf)
                
        return self.be.make_loop(self.nfpts, comm_flux)

#-------------------------------------------------------------------------------#    
    def _make_grad_at_face(self, nele):
        nvars, ndims = self.nvars, self.ndims
        lt, le, lf = self._lidx

        # Mangitude and direction of the connecting vector
        inv_tf = self._rcp_dx
        tf = self._dx_adj * inv_tf
        avec = self._vec_snorm/np.einsum('ij,ij->j', tf, self._vec_snorm)

        # Stack-allocated array
        array = self.be.local_array()

        def grad_at(i_begin, i_end, gradf, *uf):
            # Parse element views (fpts, grad)
            du = uf[:nele]
            gradu = uf[nele:]
            #print("gradf bc", gradf)
            #print("uf bc", uf)
            #print("du bc", du)
            #print("gradu bc", gradu)
            #print("le bc", le)
            #print("avec BC", avec)
            #print("tf bc", tf)
            #print("inv_tf bc", inv_tf)
            for idx in range(i_begin, i_end):
                lti, lfi, lei = lt[idx], lf[idx], le[idx]
                #*************************#
                # Complete function
                for dimension in range(ndims):
                    for variable in range(nvars):
                        gradf[dimension][variable][idx] = gradu[lti][dimension][variable][lei]
                        #gradf[dimension][variable][idx] = avec[dimension][idx] * gradu[lti][dimension][variable][lei]
                #*************************#

        return self.be.make_loop(self.nfpts, grad_at)

#-------------------------------------------------------------------------------#    
    def _make_delu(self):
        nvars = self.nvars
        lt, le, lf = self._lidx
        nf = self._vec_snorm

        bc = self.bc
        array = self.be.local_array()

        def compute_delu(i_begin, i_end, *uf):
            for idx in range(i_begin, i_end):
                ur = array(nvars)
                nfi = nf[:, idx]

                lti, lfi, lei = lt[idx], lf[idx], le[idx]
                ul = uf[lti][lfi, :, lei]
                bc(ul, ur, nfi)
                for jdx in range(nvars):
                    du = ur[jdx] - ul[jdx]
                    uf[lti][lfi, jdx, lei] = du

        return self.be.make_loop(self.nfpts, compute_delu)


#-------------------------------------------------------------------------------#    
class ParabolicDrichletBCInters(ParabolicBCInters):
    name = 'drichlet'
    def __init__(self, be, cfg, elemap, lhs, bctype):
        super().__init__(be, cfg, elemap, lhs, bctype)
        self._reqs = self.primevars


#-------------------------------------------------------------------------------#    
class ParabolicNeumannBCInters(ParabolicBCInters):
    name = 'neumann'
    def __init__(self, be, cfg, elemap, lhs, bctype):
        super().__init__(be, cfg, elemap, lhs, bctype)
        self._reqs = self.primevars
