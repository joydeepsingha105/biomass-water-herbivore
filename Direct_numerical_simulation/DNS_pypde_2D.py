# This code implements the biomass-soil_water herbivore model in two dimension
# visit https://py-pde.readthedocs.io/en/latest/getting_started.html for instructions to install py-pde and requirements with additional dependancies. 

## This is mainly used to get the Figures 12

import numba as nb  #numba
import numpy as np  #numpy
import h5py         #h5py
import matplotlib.pyplot as plt #pyplot(not used here)
from scipy.ndimage import gaussian_filter #filter to smooth out random noise (not used here) 

# all of the py-pde components used for the code. "UnitGrid" and "plot_kymograph(s)" are not used but given in case the user wish to try them out
# Please visit the py-pde documentation page for usage and function of these components. 

from pde import PDE, PDEBase, FieldCollection, UnitGrid, CartesianGrid, plot_kymograph, ScalarField, plot_kymographs, MemoryStorage, FileStorage
from pde.solvers import Controller, ExplicitMPISolver
from pde.tools import mpi
from pde.fields import VectorField
#from pde.tools.expressions import evaluate

class BWH(PDEBase):
    def __init__(self, P = 120, bc = "auto_periodic_neumann"):
        super().__init__
        self.lamda = 0.5            #vegetation growth rate
        self.E = 10.0               #root to shoot ratio
        self.K = 0.9                #maximum vegetation per unit area
        self.M = 11.4               #vegetation mortality rate
        self.DB = 1.2               #vegetation seed dispersal

        self.N = 20.0               #evaporation rate
        self.R = 0.01               #reduced evaporation due to vegetation
        self.GamW = 10.0            #soil water uptake rate
        self.DW = 150.0             #lateral soil water diffusion coefficient
        
        self.mH = 0.06              #herbivore mortality rate
        self.GamH = 0.3             #herbivore prolifereation rate
        self.AlphaH = 0.608         #vegetation consumption rate
        self.BetaH = 0.82           #satiation biomass
        self.K_H = 175.0            #maximum herbivore per unit area
        
        self.DHH = 400              #diff. coeff. for random motility
        self.Zeta = 0.001           #reference biomass for 50% random motility drop
        self.DHB = 700              #diff. coeff. for vegetaxis motility
        self.k = 0.0001             #reference biomass for 50% vegetaxis motitity drop
        
        self.P = P                  #precipitation
        self.bc = bc                #boundary condition
    
    def get_initial_state(self, grid):	 

        # Random initial condition
        B = ScalarField.random_uniform(grid)
        W = ScalarField.random_uniform(grid)
        H = ScalarField.random_uniform(grid)
        
        # changing initial condition
        B.data = 0.5 + (0.1*B.data)
        W.data = 0.8 
        H.data = 0.45
        # #this labels are useful when output data file is analysed later
        B.label = "Plants"
        W.label = "Water"
        H.label = "Herbivore"
 
        return FieldCollection([B, W, H])

    def evolution_rate(self, state, t = 0): 
        B, W, H = state
        
        #rate equations implemented in two dimension
        B_t = (self.lamda * B * W * ((1.0 + (self.E * B))**2) * (1.0 - (B/self.K))) - (self.M * B) - (self.AlphaH * H * B/(self.BetaH + B)) + (self.DB * B.laplace(self.bc))
        W_t = self.P - (self.N * W/(1.0 + (self.R * B/self.K))) - (self.GamW * B * W * (1.0 + (self.E * B))**2) + (self.DW * W.laplace(self.bc))
        H_t = (-self.mH * H) + (self.GamH * self.AlphaH * H * B * (1.0 - (H/self.K_H))/(self.BetaH + B)) - ((self.DHB * self.k/(self.k + B)) * ((H.gradient(self.bc) @ B.gradient(self.bc)))) + ((self.DHB * H * self.k/(self.k + B)**2) * (B.gradient(self.bc) @ B.gradient(self.bc))) - ((self.DHB * H * self.k/(self.k + B)) * B.laplace(self.bc)) + ((self.DHH * self.Zeta**2/(self.Zeta**2 + B**2)) * H.laplace(self.bc)) - ((self.DHH * 2 * self.Zeta**2 * B/(self.Zeta**2 + B**2)**2) * ((H.gradient(self.bc) @ B.gradient(self.bc)))) 
        return FieldCollection([B_t, W_t, H_t])
    

    #numba implemention to run this code parallely. Python packages such as Numba, Numba-mpi, mpi4py are required to
    #run it parallely. We suggest to consult the py-pde documentation page mentioned in the beginning.
    def _make_pde_rhs_numba(self, state):
        """numba implementation of the PDE"""
        lamda = self.lamda
        E = self.E
        K = self.K
        M = self.M
        DB = self.DB

        N = self.N
        R = self.R
        GamW = self.GamW
        DW = self.DW
        
        mH = self.mH
        GamH = self.GamH
        AlphaH = self.AlphaH
        BetaH = self.BetaH
        K_H = self.K_H
        
        DHH = self.DHH
        Zeta = self.Zeta
        DHB = self.DHB
        k = self.k
        
        P = self.P

        #defining the gradient, square of gradient and laplace operation
        #these are defined according to py-pde methods. Please consult the py-pde documentation to get familer.
        gradient = state.grid.make_operator("gradient", bc = self.bc)
        gradient_squared = state.grid.make_operator("gradient_squared", bc = self.bc)
        laplace = state.grid.make_operator("laplace", bc = self.bc) 
        dot = VectorField(state.grid).make_dot_operator()
        @nb.njit
        def pde_rhs(state_data, t):
            B = state_data[0]
            W = state_data[1]
            H = state_data[2]
            
            rate = np.empty_like(state_data)
            #rate equations implemented in two dimension
            rate[0] = (lamda * B * W * ((1.0 + (E * B))**2) * (1.0 - (B/K))) - (M * B) - (AlphaH * H * B/(BetaH + B)) + (DB * laplace(B))
            rate[1] = P - (N * W/(1.0 + (R * B/K))) - (GamW * B * W * (1.0 + (E * B))**2) + (DW * laplace(W))
            rate[2] = (-mH * H) + (GamH * AlphaH * H * B * (1 - (H/K_H))/(BetaH + B)) - ((DHB * k/(k + B)) * (dot(gradient(H), gradient(B)))) + ((DHB * H * k/(k + B)**2) * gradient_squared(B)) - ((DHB * H * k/(k + B)) * laplace(B)) + (((DHH * Zeta**2)/((Zeta**2) + (B**2))) * laplace(H)) - ((DHH * 2 * (Zeta**2) * B/((Zeta**2) + (B**2))**2) * dot(gradient(H), gradient(B)))
            return rate
        
        return pde_rhs

#domain size, same domain size is used in both x and y directions 
Lx = 16*np.pi

#mesh size, same mesh size is used in both x and y directions
mesh = 512

#defining the grid
grid = CartesianGrid([[0,Lx],[0,Lx]], [mesh,mesh], periodic = [True,True])



P_list = [345]          #precipitation
#useful when using jupyter notebook for the analysis of data without saving a data file
storage = MemoryStorage()

for p in P_list:
    eq = BWH()      #initialising the PDE class
    eq.P = p        #P input 
    storage = FileStorage("2D_data.h5")     #output data file
    trackers = ["progress", storage.tracker(1)]     #tracker to see the duration of simulation 

    #defining the solver : Euler  with adaptive time stepping with domain divided into 32 cores
    solver_mpi = ExplicitMPISolver(eq, scheme = 'euler', decomposition=32, backend = 'auto', adaptive = True)
    
    #controller for parallel computation specifying the total time duration 
    controller_mpi = Controller(solver_mpi, t_range = 2000, tracker = trackers)

    #initial condition
    state = eq.get_initial_state(grid);	

    #running of the model 
    sol = controller_mpi.run(state)
