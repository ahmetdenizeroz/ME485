# -*- coding: utf-8 -*-
from solvers.base.system import BaseSystem
from solvers.advection import AdvectionElements, AdvectionIntInters, AdvectionMPIInters,AdvectionBCInters, AdvectionVertex
import sys

class AdvectionSystem(BaseSystem):
    name = 'advection'
    _elements_cls = AdvectionElements
    _intinters_cls = AdvectionIntInters
    _bcinters_cls = AdvectionBCInters
    _mpiinters_cls = AdvectionMPIInters
    _vertex_cls = AdvectionVertex

    def rhside(self, idx_in=0, idx_out=1, t=0, is_norm=True):
         # Adjust Banks
        self.eles.upts_in.idx = idx_in
        self.eles.upts_out.idx = idx_out

        # Queue for MPI
        q = self._queue

        # Compute solution at flux point (face center)
        self.eles.compute_fpts()
        print("hear1")
        if self.mpiint:
            # Start MPI communication for Inters
            self.mpiint.pack()
            self.mpiint.send(q)
            self.mpiint.recv(q)
        print("hear2")
        # Compute Difference of solution at Inters
        self.iint.compute_delu()
        print("hear3")
        self.bint.compute_delu()
        print("hear4")
        if self.mpiint:
            # Finalize MPI communication
            q.sync()

            # Compute Difference of solution at MPI Inters
            self.mpiint.compute_delu()
        
        self.eles.compute_grad()
        print("hear5")
        if self._limiter=='mlp-u1' or self._limiter=='mlp-u1':
            # Compute extreme values at vertex
            self.vertex.compute_extv()
            if self.vertex.mpi:
                # Start MPI communication for Vertex
                self.vertex.pack()
                self.vertex.send(q)
                self.vertex.recv(q)
            if self.vertex.mpi:
                # Finalize MPI communication
                q.sync()
                # Unpack (Sort vetex extremes)
                self.vertex.unpack()
        elif self._limiter=='barth-jespersen' or  self._limiter=='venkatrishnan' :
            # Compute solution at flux point (face center)
            self.eles.compute_fpts()
            self.iint.compute_minmax()
            self.bint.compute_minmax()

        print("hear6")
        # Compute slope limiter
        self.eles.compute_limiter()
        print("hear7")
        # print(self.eles.lim)
        # Compute reconstruction
        self.eles.compute_recon()
        print("hear8")
        if self._is_recon and self.mpiint:
            # Start MPI communication to exchange reconstructed values at face
            self.mpiint.pack()
            self.mpiint.send(q)
            self.mpiint.recv(q)
        print("hear9")
        # # Compute flux
        self.iint.compute_flux()
        print("hear9.5")
        self.bint.compute_flux()
        print("hear10")
        if self.mpiint:
            # Finalize MPI communication
            q.sync()

            # Compute flux at MPI Inters
            self.mpiint.compute_flux()
        print("hear11")
        # Compute divergence 
        self.eles.div_upts(t)
        print("hear12")
        if is_norm:
            # Compute residual if requested
            resid = sum(self.eles.compute_norm())
            return resid
        else:
            return 'none'

        # sys.exit()
        print("hear13")
#-------------------------------------------------------------------------------#    
    def timestep(self, cfl, idx_in=0):
        # Compute time step with the given CFL number
        self.eles.upts_in.idx = idx_in
        self.eles.timestep(cfl)

#-------------------------------------------------------------------------------#
    def post(self, idx_in=0):
        # Post-process
        self.eles.upts_in.idx = idx_in
        self.eles.post()
