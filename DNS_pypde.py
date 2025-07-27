import numba as nb 
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from pde import PDE, PDEBase, FieldCollection, UnitGrid, CartesianGrid, plot_kymograph, ScalarField, plot_kymographs, MemoryStorage, FileStorage
from pde.solvers import Controller, ExplicitMPISolver
from pde.tools import mpi 
class BWH(PDEBase):
    def __init__(self, P = 120, bc = "auto_periodic_neumann"):
        super().__init__
        self.lamda = 0.5            #vegetation growth rate
        self.E = 10.0               #root to shoot ratio
        self.K = 0.9                #maximal vegetation capacity per unit square
        self.M = 11.4               #vegetation mortality rate
        self.DB = 1.2               #vegetation diffusion representing seed dispersal

        self.N = 20.0               #evaporation rate
        self.R = 0.01               #reduced evaporation due to shading
        self.GamW = 10.0            #soil water uptake rate
        self.DW = 150.0             #lateral diffusion coefficient of soil water
        
        self.mH = 0.06              #mortality rate of herbivores
        self.GamH = 0.3             #proliferation rate of herbivores 
        self.AlphaH = 0.608         #consumption rate of herbivores
        self.BetaH = 0.82           #satiation biomass
        
        self.DHH = 400              #maximal random motility of herbivores
        self.Zeta = 0.001           #reference vegetation-biomass density for 50% motility drop
        self.DHB = 700              #maximal vegetaxis motility 
        self.k = 0.0001             #reference vegetation-biomass density for 50% motility drop
        
        self.P = P                  #precipitation rate
        self.bc = bc                #boundary condition
    
    def get_initial_state(self, grid):                  #this function initialises the variables
        # random initial condition
        B = ScalarField.random_uniform(grid)            #vegetation 
        W = ScalarField.random_uniform(grid)            #soil water 
        H = ScalarField.random_uniform(grid)            #herbivores 
        #change the initial condition and update in B.data, W.data, H.data

        B.label = "Plants"                              #labelling the variables for the data file 
        W.label = "Water"
        H.label = "Herbivore"
 
        return FieldCollection([B, W, H])

    def evolution_rate(self, state, t = 0):             #description of the rate equations
        B, W, H = state
        
        #vegetation
        B_t = (self.lamda * B * W * ((1.0 + (self.E * B))**2) * (1.0 - (B/self.K))) - (self.M * B) - (self.AlphaH * H * B/(self.BetaH + B)) + (self.DB * B.laplace(self.bc))
        
        #soil water
        W_t = self.P - (self.N * W/(1.0 + (self.R * B/self.K))) - (self.GamW * B * W * (1.0 + (self.E * B))**2) + (self.DW * W.laplace(self.bc))
        
        #herbivore
        H_t = (-self.mH * H) + (self.GamH * self.AlphaH * H * B/(self.BetaH + B)) - ((self.DHB * self.k/(self.k + B)) * ((H.gradient(self.bc) * B.gradient(self.bc)))) + ((self.DHB * H * self.k/(self.k + B)**2) * B.gradient(self.bc)**2) - ((self.DHB * H * self.k/(self.k + B)) * B.laplace(self.bc)) + ((self.DHH * self.Zeta**2/(self.Zeta**2 + B**2)) * H.laplace(self.bc)) - ((self.DHH * 2 * self.Zeta**2 * B/(self.Zeta**2 + B**2)**2) * ((H.gradient(self.bc) * B.gradient(self.bc))))
            
        return FieldCollection([B_t, W_t, H_t])
    
    def _make_pde_rhs_numba(self, state):   #numba implementation for parallelisation
        lamda = self.lamda                  #growth rate of vegetation  
        E = self.E                          #root to shoot ratio
        K = self.K                          #maximal vegetation capacity per unit square
        M = self.M                          #mortality rate of vegetation
        DB = self.DB                        #vegetation diffusion coefficient due to seed dispersal

        N = self.N                          #soil water evaporation rate
        R = self.R                          #reduced evaporation due to vegetation
        GamW = self.GamW                    #soil wate uptake rate 
        DW = self.DW                        #lateral diffusion coefficient of soil water
        
        mH = self.mH                        #mortality rate of herbivores
        GamH = self.GamH                    #proliferation rate of herbivores 
        AlphaH = self.AlphaH                #consumption rate of herbivores
        BetaH = self.BetaH                  #satiation biomass
        
        DHH = self.DHH                      #maximal random motility
        Zeta = self.Zeta                    #reference vegetation biomass density for 50% motility drop
        DHB = self.DHB                      #maximal vegetaxis motility
        k = self.k                          #reference vegetation biomass density for 50% motility drop
        
        P = self.P                          #precipitation
        
        #defining the gradient, square of gradient and the laplace operation with periodic boundary condition
        gradient = state.grid.make_operator("gradient", bc = self.bc)                   
        gradient_squared = state.grid.make_operator("gradient_squared", bc = self.bc)
        laplace = state.grid.make_operator("laplace", bc = self.bc)
        
        @nb.jit
        def pde_rhs(state_data, t):         #function for rate equations of the variables
            B = state_data[0]               #vegetation
            W = state_data[1]               #soil water
            H = state_data[2]               #herbivore
            
            rate = np.empty_like(state_data)

            #vegetation
            rate[0] = (lamda * B * W * ((1.0 + (E * B))**2) * (1.0 - (B/K))) - (M * B) - (AlphaH * H * B/(BetaH + B)) + (DB * laplace(B))
            
            #soil water
            rate[1] = P - (N * W/(1.0 + (R * B/K))) - (GamW * B * W * (1.0 + (E * B))**2) + (DW * laplace(W))
            
            #herbivore
            rate[2] = (-mH * H) + (GamH * AlphaH * H * B/(BetaH + B)) - ((DHB * k/(k + B)) * ((gradient(H) * gradient(B)))) + ((DHB * H * k/(k + B)**2) * gradient_squared(B)) - ((DHB * H * k/(k + B)) * laplace(B)) + (((DHH * Zeta**2)/((Zeta**2) + (B**2))) * laplace(H)) - ((DHH * 2 * (Zeta**2) * B/((Zeta**2) + (B**2))**2) * (gradient(H) * gradient(B)))
            return rate
        
        return pde_rhs


Lx = 2*np.pi                                                #domain length
mesh = 256                                                  #grid size
grid = CartesianGrid([[0,Lx]], [mesh], periodic = True)     #setting grid type and periodic boundary condition

P_list = [250]                                              #setting the precipitation value                                   
for p in P_list:                                            #calling the BWH class for each precipitation value 
    eq = BWH()                                              #initialting the BWH class
    eq.P = p                                                #precipitation
    storage = FileStorage("TW_H.h5")                        #output data file
    trackers = ["progress", storage.tracker(0.1)]           #this will track how much time is left during the run
    solver_mpi = ExplicitMPISolver(eq, scheme = 'rk', decomposition=2, backend = 'auto', adaptive = True)   #Mpi solver to specify the adaptive time stepper, core usage
    controller_mpi = Controller(solver_mpi, t_range = 2000, tracker = trackers)                             #total time duration of the run 
    state = eq.get_initial_state(grid)	                                                                    #initialising the grid and initial condition 
    sol = controller_mpi.run(state)                                                                         #final run command


