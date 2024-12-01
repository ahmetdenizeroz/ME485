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
    def compute_L2_norm(self, exact_soln):
        nvars, nface, ndims = self.nvars, self.nface, self.ndims
        vol                 = self._vol
         # **************************#

        neles = self.neles

        resid = 0

        #Assume that there is a matrix of results. The dimensions of this matrix are num.of.elements X num.of.variables . 
        #Example: Gradient vectors at the center of elements. (1st row is for element1, 2nd row is for element2, etc.)
        #The following code calculates the weighted L2 error norm of this matrix of results.

        for i in range(neles):  #Loop over elements (rows of the matrix)
            for j in range(nvars):  #Loop over variables (columns of the matrix)
                
                #Finite volume method L2 formula. This formula takes squares of variables which belong to the same element, sum them and
                #multiply the result with the volume of the element. Therfore, the sum of square of variables are weighted wrt. the size of
                #the element. Then, all the weighted sum of squares are again summed up for different elements and stored in the resid value.
                #The summations value is divided to the size of matrix to find an average value. Finally, the sqrt of the value is taken.

                resid = resid + np.sum( ( self.upts_out[j, i] - exact_soln[j,i] ) ** 2 * vol[i] )  

        resid = resid / (neles * nvars)
        resid = np.sqrt(resid)
    
        # **************************#

        return resid

#-------------------------------------------------------------------------------#
    # Assign cell centers values to face centers
    def _make_compute_fpts(self):
        nvars, nface = self.nvars, self.nface

        def _compute_fpts(i_begin, i_end, upts, fpts):
            # Copy upts to fpts
            # fpts --> Flux Points
            # upts --> Universal Points (Center Points)

            for idx in range(i_begin, i_end): #Iteration over elements
                for j in range(nvars):  #Iteration over variables
                    # **************************#
                    for k in range(nface):   #Iteration over faces
                        
                        #The values calculated for the center (upts) are stored in the faces of the element.
                        fpts[k, j, idx] = upts[j, idx]    

                    # **************************#
        
        return self.be.make_loop(self.neles, _compute_fpts)   
#-------------------------------------------------------------------------------#
    def _make_grad_ls(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        # Gradient operator 
        op = self._grad_operator   #[face, element]

        dxf = self.dxf  #something new is added here! [faces, elem, dim](?)
        
        #Assume that the matrix system is A*x = b
        
        def _cal_grad(i_begin, i_end, fpts, grad):
            # Elementwiseloop
            
            grad = []
            for element_index in range(i_begin, i_end):  #Loop over elements
                # **************************#
                grad_forelement = []
                for variable_index in range(nvars):  #Loop over variables

                    A = []
                    b = []
                    for face_index in range(nface):   #Loop over dimensions

                        A_oneface = []
                        for dim_index in range(ndims):    #Loop over faces

                            A_oneface = A_oneface + [dxf(face_index, element_index, dim_index)] 
                            #dxf gives the value of self.xf - self.xc (position vector between face and centroid)
                            #A_oneface creates an array of [x1_x0 y1-y0 ....] then [x2_x0 y2-y0 ....] ...

                    b = np.append(b, op(face_index, element_index) * [ fpts[face_index, variable_index, element_index] - self.upts_out[variable_index, element_index] ], axis = 0)
                    # The b vector which is [w1*[Q1 - Qc], w2*[Q2 - Qc], ...] (w --> weight factor which came from op() )

                    A = np.append(A, A_oneface, axis = 0) 
                    #The arrays are combined and a matrix is created. 
                    #Rows are faces, columns are dimensional pos. differences.

                grad_forelement = np.append(grad_forelement, np.linalg.inv(np.transpose(A) * A ) * np.transpose(A) * b, axis = 0)
            
            grad = np.append(grad, grad_forelement, axis = 2)

                # **************************#

        # Compile the function
        return self.be.make_loop(self.neles, _cal_grad)    
#-------------------------------------------------------------------------------#
    def _make_grad_gg(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        #Normal vector and volume
        snorm_mag = self._mag_snorm_fpts
        snorm_vec = np.rollaxis(self._vec_snorm_fpts, 2)
        vol       = self._vol

        def _cal_grad(i_begin, i_end, fpts, grad):
            # Elementwise loop starts here
            for i in range(i_begin, i_end):
               
                # **************************#
                
                for j in range(nvars):  

                    for k in range(ndims):  

                        grad[k, j, i] = 0.0 

                        for m in range(nface):  

                            grad[k, j, i] = grad[k, j, i] + snorm_vec[k, m, i] * fpts[m, j, i] / vol[i]


                # **************************#

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
        grad_method = self._grad_method
        
        op = []

        p = 1
        
      
        if grad_method == 'weighted-least-square':

            for element_index in range(self.nelem):

                for face_index in range(self.nface):

                    op_element = np.append(op, 1/(distance(1, face_index, element_index))**p , axis = 0)
            
            op = op + op_element  # op => [face, element]

        else:
            
            op = np.ones_like(distance)

        # **************************#
        return op

    def _make_post(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        grad = self.grad
        xc = self.xc.T        
        def post(i_begin, i_end, upts):
            # Apply the function over elements
            for idx in range(i_begin, i_end):          
                print(idx, grad[0, 0, idx], grad[1, 0, idx])

        return self.be.make_loop(self.neles, post)
