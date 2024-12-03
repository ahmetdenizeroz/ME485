# -*- coding: utf-8 -*-
import numpy as np
import re

from solvers.base import BaseElements, BaseIntInters
from backends.types import ArrayBank, Kernel, NullKernel
from utils.np import eps, chop, npeval
from solvers.grad.inters import GradIntInters

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
        
        # Assuming lhs and rhs are determined somewhere or passed to the constructor
        # You need to pass the correct `lhs` and `rhs` for the face connectivity
        lhs = self._generate_lhs(eles)  # Method to generate lhs (list of element-face-point tuples)
        rhs = self._generate_rhs(eles)  # Method to generate rhs (list of element-face-point tuples)
        
        # Create an instance of GradIntInters, passing lhs and rhs
        self.grad_int_inters = GradIntInters(be, cfg, eles, lhs, rhs)  # Initialize with lhs and rhs
        
    def _generate_lhs(self, eles):
        # Generate the lhs list of element-face-point tuples based on the provided elements
        lhs = []
        for ele in eles:
            # Create a tuple for each element, face, and point combination (example)
            for face in ele.faces:
                for point in face.points:
                    lhs.append((ele.type, ele.id, face.id, point.id))
        return lhs

    def _generate_rhs(self, eles):
        # Generate the rhs list of element-face-point tuples, which should be computed
        # based on the neighboring element's faces or some other criteria
        rhs = []
        for ele in eles:
            # Assuming rhs is similar to lhs but involves neighboring elements or different logic
            for face in ele.neighbors:
                for point in face.points:
                    rhs.append((ele.type, ele.id, face.id, point.id))
        return rhs        
        
    
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
        upts = self.upts_in  # Computed solution at each element

        # Define the exact solution function
        def exact_solution(x):
            # Example: quadratic solution x^2 (for a single variable case)
            return x**2

        # Compute the exact solution at the element centers
        exact_upts = np.zeros_like(upts)
        for i in range(upts.shape[-1]):  # Loop over elements
            # Assuming `xc` contains the coordinates of element centers
            x = self.xc[:, i]  # Element center coordinates
            exact_upts[:, i] = exact_solution(x)

        # Compute the squared error
        squared_error = np.zeros(upts.shape[-1])  # One error per element
        for i in range(upts.shape[-1]):  # Loop over elements
            for l in range(nvars):  # Loop over variables
                squared_error[i] += np.sum((exact_upts[l, i] - upts[l, i])**2)

        # Weight by volume and compute the total L2 error
        weighted_error = squared_error * vol
        total_error = np.sum(weighted_error)

        # Normalize by total volume
        total_volume = np.sum(vol)
        resid = np.sqrt(total_error / total_volume)

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
                for l in range(nvars):  # Loop over variables
                    # Extract the operator (A matrix) for the current element
                    A = op[:, :, i]  # Shape: (Nfaces, Ndims)

                    # Construct b vector for the current variable and element
                    # b = phi_f - phi_center
                    phi_f = fpts[:, l, i]  # Shape: (Nfaces,)
                    phi_center = self.upts[l, i]  # Use `upts` for the central value
                    b = phi_f - phi_center  # Shape: (Nfaces,)

                    # Solve for the gradient: grad = (A^T A)^-1 A^T b
                    # Using np.linalg.lstsq for least-squares solution
                    grad[:, l, i], *_ = np.linalg.lstsq(A, b, rcond=None)
        return self.be.make_loop(self.neles, _cal_grad)       
#-------------------------------------------------------------------------------#
    def _make_grad_gg(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        # Normal vector and volume
        snorm_mag = self._mag_snorm_fpts
        snorm_vec = np.rollaxis(self._vec_snorm_fpts, 2)
        vol       = self._vol
        # Precompute face averages and differences
        avgu = self.grad_int_inters.compute_avgu
        delu = self.grad_int_inters.compute_delu
     
        
        def _cal_grad(i_begin, i_end, fpts, grad):
            for i in range(i_begin, i_end):  # Loop through elements
                for l in range(nvars):  # Loop through variables
                    grad_val = np.zeros(ndims)  # Initialize gradient vector
                
                    # Compute face-centered values using midpoint formula
                    for f in range(nface):
                        # Fetch precomputed averaged face value (phi_f) Note that w = 0.5
                        phi_f = avgu[f, l, i]

                        # Accumulate the gradient contributions
                        for d in range(ndims):
                            grad_val[d] += phi_f * snorm_vec[f, d, i] * snorm_mag[f, i]

                    # Normalize by volume to get the gradient
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
    
        # Get necessary attributes
        ndims = self.ndims       # Number of spatial dimensions
        nface = self.nface       # Number of faces per element
        neles = self.neles       # Total number of elements
        phi = self.upts_in       # Primitive variable values at element centers
        dxc = np.rollaxis(self.dxc, 2)  # Rearrange axes (Nelements, Nfaces, Ndims)

        # Precompute distances
        distance = np.linalg.norm(dxc, axis=-1)  # (Nelements, Nfaces) !!! Orjinalinde axis = 0 vermiş 

        # Create gradient operator matrix
        op = np.zeros((ndims, nface, neles))  # Placeholder

        # Loop over elements
        for e0 in range(neles):
            # Collect neighbors' displacement vectors and distances
            rel_pos = dxc[e0]  # Shape: (Nfaces, Ndims)
            dist = distances[e0]  # Shape: (Nfaces,)

            # Collect variable differences
            phi_e0 = phi[e0]
            phi_neighbors = self.get_neighbors_phi(e0)  # Fetch variable values at neighbors
            b = phi_neighbors - phi_e0  # Shape: (Nfaces,)

            # Construct weight matrix if weighted LS is used
            p = 1  # Default p=0
            if p > 0:
                W = np.diag(1.0 / dist**p)
                A = rel_pos
                x = np.linalg.inv(A.T @ W @ A) @ (A.T @ W @ b)
            else:
                A = rel_pos
                x = np.linalg.inv(A.T @ A) @ (A.T @ b)

            # Store computed gradient
            op[:, :, e0] = x.reshape((ndims, -1))

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
