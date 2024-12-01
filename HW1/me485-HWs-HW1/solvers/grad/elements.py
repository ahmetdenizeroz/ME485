    # -*- coding: utf-8 -*-
import numpy as np
import re

from solvers.base import BaseElements
from backends.types import ArrayBank, Kernel, NullKernel
from utils.np import eps, chop, npeval


class gradFluidElements:
    @property
    def primevars(self):
        # Primitive variables
        return ['q'] 
    @property
    def conservars(self):
        # Conservative variables
        pri = self.primevars

        # rho,rhou,rhov,rhow,E
        return [pri[0]]

    def prim_to_conv(self, pri, cfg):
        return [pri[0]]

    def conv_to_prim(self, con, cfg):
        return [con[0]]

    def fix_nonPys_container(self):
        # Constants and dimensions
        ndims, nfvars = self.ndims, self.nfvars

        def fix_nonPhy(u):
            u[0] = u[0]
        # Compile the function
        return self.be.compile(fix_nonPhy)

class GradElements(BaseElements,  gradFluidElements):
    nreg = 1
    def __init__(self, be, cfg, name, eles):
        super().__init__(be, cfg, name, eles)
        self.nvars = len(self.primevars)
        #print(self.nvars)
        self.nfvars = self.nvars
        self._const = cfg.items('constants')
        self._grad_method = cfg.get('solver', 'gradient')
        # print(self.dxv)

    def construct_kernels(self, vertex, nreg, impl_op=0):
        self.vertex = vertex

        # Upts : Solution vector
        self.upts = upts = [self._ics.copy() for i in range(nreg)]
        del(self._ics)

        # Solution vector bank and assign upts index
        self.upts_in = upts_in = ArrayBank(upts, 0)
        self.upts_out = upts_out = ArrayBank(upts, 1)
        #print(upts_in)
        #print(upts_out)

        # Construct arrays for flux points, dt and derivatives of source term
        self.fpts = fpts = np.empty((self.nface, self.nvars, self.neles))
        self.grad = grad = np.zeros((self.ndims, self.nvars, self.neles))
        print(np.shape(fpts))
        print(np.shape(fpts))

        lim = np.ones((self.nvars, self.neles))

        limiter = self.cfg.get('solver', 'limiter', 'none')
        # Prepare vertex array
        vpts = vertex.make_array(limiter)
        
        self.compute_fpts = Kernel(self._make_compute_fpts(), upts_in, fpts)

       
        if( (self._grad_method=='least-square') or (self._grad_method=='weighted-least-square')):
            self.compute_grad = Kernel(self._make_grad_ls(), fpts, grad)       
        elif(self._grad_method == 'green-gauss-cell'):
            self.compute_grad = Kernel(self._make_grad_gg(), fpts, grad)
        else:
            self.compute_grad = Kernel(self._make_grad_ls(), fpts, grad)

        
        # Kernel for linear reconstruction
        self.compute_recon = Kernel(self._make_recon(), upts_in, grad, lim, fpts)
        
        # Kernel to post-process
        self.post = Kernel(self._make_post(), upts_in)

 #-------------------------------------------------------------------------------#
    def compute_L2_norm(self):
        nvars, nface, ndims = self.nvars, self.nface, self.nvars
        vol                 = self._vol

        # **************************#

        # **************************#

        return resid

#-------------------------------------------------------------------------------#
    # Assign cell centers values to face centers
    def _make_compute_fpts(self):
        nvars, nface = self.nvars, self.nface

        # UPTS is an array containing center values of every element in the region of interest.
        # FPTS is an array containing face values. It's sahepe is (nface, nvars, nelemets), menaing that
        # it can store face values of every variable for every element, for numbe of face times. Meaning that
        # different values of coincident faces can ve stored for different elements.

        def _compute_fpts(i_begin, i_end, upts, fpts):
            # Copy upts to fpts
            for idx in range(i_begin, i_end):
                for j in range(nvars):
                    for face in range(nface):
                        fpts[face, 0, idx] = upts[0, 0, idx]
        return self.be.make_loop(self.neles, _compute_fpts)

#-------------------------------------------------------------------------------#
    def _make_grad_ls(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        # Gradient operator 
        op = self._grad_operator

        def _cal_grad(i_begin, i_end, fpts, grad):
            # Elementwiseloop
            for i in range(i_begin, i_end):
                for variable in range(nvars):
                    b = np.zeros((nface, 1))
                    for face in range(nface):
                        b[face] = fpts[face, variable, i]
                    grad = op[i] @ b
                    for dimension in range(ndims):
                        grad[dimension, variable, i] = grad[dimension]
        # Compile the function
        return self.be.make_loop(self.neles, _cal_grad)

#-------------------------------------------------------------------------------#
    def _make_grad_gg(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        #Normal vector and volume
        snorm_mag = self._mag_snorm_fpts
        snorm_vec = np.rollaxis(self._vec_snorm_fpts, 2)
        vol       = self._vol

        # The main consideretion of impelentation of green gauss cell based method is that how to reach face values of
        # concidend faces. After spending some time in the solver, i decided that _make_avgu methods finds the avarage
        # values on the faces, and updates the fpts accordingly.

        def _cal_grad(i_begin, i_end, fpts, grad):
            # Elementwise loop starts here
            for i in range(i_begin, i_end):
                for dimension in range(ndims):
                    RHS = 0
                    for face in range(nface):
                        RHS += fpts[face, 0, i] * snorm_vec[dimension, face, i] * snorm_mag[face, i]
                    grad[dimension, 0, i] = (1/vol) * RHS

        # Compile the function
        return self.be.make_loop(self.neles, _cal_grad)      

#-------------------------------------------------------------------------------#
    def _make_recon(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars

        # Displacement vector
        op = self.dxf

        def _cal_recon(i_begin, i_end, upts, grad, lim, fpts):
            # Elementwise dot product and scale with limiter
            # TODO: Reduce accesing global array
            for i in range(i_begin, i_end):
                for l in range(nvars):
                    for k in range(nface):
                        tmp = 0
                        for j in range(ndims):
                            tmp += op[k, j, i] * grad[j, l, i]
                                                    
                        fpts[k, l, i] = upts[l, i] + lim[l, i]*tmp

        return self.be.make_loop(self.neles, _cal_recon)


#-------------------------------------------------------------------------------#
    @property
    # @fc.lru_cache()
    # @chop
    def _grad_operator(self):
        # Difference of displacement vector (cell to cell)
        # (Nfaces, Nelements, dim) -> (dim, Nfaces, Nelements)
        dxc = np.rollaxis(self.dxc, 2)
        # (Nfaces, Nelements)

        # only possible use of distance is for weighted ls.
        distance = np.linalg.norm(dxc, axis=0)

        # we need p value assignment for determining the weight, it will be taken as 2 for now.
        # with examination of code, a little help is taken from base\elements.py
        p = 2
        w = 1

        if self._grad_method == 'least-square':
            w = 1.0
        elif self._grad_method == 'weighted-least-square':
            w = 1 / (distance**p)

        # Scaled dxc vector
        dxcs = dxc*np.sqrt(w)

        # Creating empty array for operator storage
        op = np.empty(self.neles, dtype=object)
        for element in range(self.neles):
            A = np.zeros((self.nface, self.ndims))
            for face in range(self.nface):
                for dimension in range(self.ndims):
                    A[face][dimension] = dxcs[dimension][face][element]
            AT = np.transpose(A)
            mult = np.dot(AT , A)
            op[element] = np.linalg.inv(mult) @ AT

        return op

    def _make_post(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        grad = self.grad
        xc = self.xc.T        
        def post(i_begin, i_end, upts):
            # Apply the function over eleemnts
            for idx in range(i_begin, i_end):          
                print(idx, grad[0, 0, idx], grad[1, 0, idx])

        return self.be.make_loop(self.neles, post)