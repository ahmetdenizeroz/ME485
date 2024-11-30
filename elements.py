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

        # Construct arrays for flux points, dt and derivatives of source term
        self.fpts = fpts = np.empty((self.nface, self.nvars, self.neles))
        self.grad = grad = np.zeros((self.ndims, self.nvars, self.neles))

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
        nvars, nface, ndims = self.nvars, self.nface, self.ndims
        vol                 = self._vol
        # **************************#

        #Complete

        # **************************#

        return resid

#-------------------------------------------------------------------------------#
    # Assign cell centers values to face centers
    def _make_compute_fpts(self):
        nvars, nface = self.nvars, self.nface

        def _compute_fpts(i_begin, i_end, upts, fpts):
            # Copy upts to fpts
            for idx in range(i_begin, i_end): # Loop through each element in the given range
                for j in range(nvars): # Loop through each variable
                    # **************************#

                    for f in range(nface): # Loop through each face
                        fpts[f, j, idx] = upts[j, idx]

                    # **************************#       
        return self.be.make_loop(self.neles, _compute_fpts)   
#-------------------------------------------------------------------------------#
    def _make_grad_ls(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        # Gradient operator 
        op = self._grad_operator
        
        def _cal_grad(i_begin, i_end, fpts, grad):
            # Elementwiseloop
            for i in range(i_begin, i_end):
                # **************************#
                #Complete
                # **************************#
                # Compile the function
        return self.be.make_loop(self.neles, _cal_grad)    
#-------------------------------------------------------------------------------#
    def _make_grad_gg(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        # Normal vector and volume
        snorm_mag = self._mag_snorm_fpts
        snorm_vec = np.rollaxis(self._vec_snorm_fpts, 2)
        vol       = self._vol
        r_centroid = self._r_centroid  # Element centroids
        r_face = self._r_face  # Face positions

        def _cal_grad(i_begin, i_end, fpts, grad):
            max_iter = 10  # Maximum iterations for convergence
            tol = 1e-6  # Tolerance for convergence

            for i in range(i_begin, i_end):  # Loop through elements
                for l in range(nvars):  # Loop through variables
                    grad_val = np.zeros(ndims)  # Initialize gradient vector
                    phi_face = np.zeros(nface)  # Initialize face values

                    # Step 1: Initial midpoint method for face values
                    for f in range(nface):
                        adj_elem = self.get_adjacent_element(i, f)  # Get adjacent element index
                        if adj_elem is not None:
                            phi_face[f] = 0.5 * (fpts[f, l, i] + fpts[f, l, adj_elem])
                        else:
                            phi_face[f] = fpts[f, l, i]  # Boundary face

                    # Step 2: Iterative correction
                    for k in range(max_iter):
                        phi_face_new = np.copy(phi_face)
                        for f in range(nface):
                            adj_elem = self.get_adjacent_element(i, f)
                            if adj_elem is not None:
                                r_e = r_centroid[:, i]
                                r_ef = r_centroid[:, adj_elem]
                                r_f = r_face[:, f]

                                grad_e = grad[:, l, i]
                                grad_ef = grad[:, l, adj_elem]

                                correction = 0.5 * (grad_e + grad_ef) @ (r_f - 0.5 * (r_e + r_ef))
                                phi_face_new[f] = phi_face[f] + correction

                        # Check for convergence
                        if np.linalg.norm(phi_face_new - phi_face) < tol:
                            break
                        phi_face = phi_face_new

                    # Step 3: Compute the Green-Gauss gradient
                    for f in range(nface):
                        for d in range(ndims):
                            grad_val[d] += phi_face[f] * snorm_vec[f, d, i] * snorm_mag[f, i]

                    # Normalize by volume
                    for d in range(ndims):
                        grad[d, l, i] = grad_val[d] / vol[i]


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
                            tmp += op[k, j, i]*grad[j, l, i]
                                                    
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
        distance = np.linalg.norm(dxc, axis=0)


        # **************************#
        #Complete
        # **************************#
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
