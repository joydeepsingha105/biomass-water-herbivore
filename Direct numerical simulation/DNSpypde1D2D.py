# This code implements the biomass-soil_water herbivore model in one dimension and comments are provided for easily implenting it two dimension. 

# visit https://py-pde.readthedocs.io/en/latest/getting_started.html for instructions to install py-pde and requirements with additional dependancies. 

# This is mainly used to get the Figures 8, 11, 13

# Figure 8 is obtaine from simulation of one dimensional implementation where the initial condition was a sum of a randon noise and constant throughout the domain. The default output file is h5py format. Plot script are in pyplot which are available upon reasonable request to the authors.  

# Figure 11 is obtained step by step for the one dimensional case. First the precipitation range is identified where a single patch of vegetation spot is stable. 
# The conventional way to do it is to first obtain a patterned solution in the "B, W" system without "H". The code can be easily changed to a simulate "B,W" model by commenting out the "H" equation and removing the consumption term from the "B" rate equations. From the output data file, a domain length of size of a wave length containing one containing one wave was chosen and kept as it is while the data for the rest of the domain was changed to zero. This modified data file was used as initial condition for "B, W" system with a slightly lower P. Thus this process was repeated until the single patch intitial condition was stable as it is. During this process we reduced P by a value of 2 in each trial step until the range isola was obtained. The single spot solution then was used as an initial condition for the "B, W, H". Initial condition for the "H" was done through a Gaussian function (shown here in the relevant section of the code).   

# Figure 13 is obtained in similar step as in Figure 11. However here we first obtain the spot pattern in the "B, W" model. In the same way as before we isolate a spot and identify P value where only this spot are stable. We then introduce "H" through a Gaussian function for the complete "B, W, H" where the stable single spot solution was used as an initial condition for B ans W.   

import numba as nb      #numba
import numpy as np      #numpy
import h5py             #h5py 
import matplotlib.pyplot as plt     #matplotlib (not used here)
from scipy.ndimage import gaussian_filter  #to smooth out random initial condition (not used here)

# all of the py-pde components used for the code. "UnitGrid" and "plot_kymograph(s)" are not used but given in case the user wish to try them out
# Please visit the py-pde documentation page for usage and function of these components. 


from pde import PDE, PDEBase, FieldCollection, UnitGrid, CartesianGrid, plot_kymograph, ScalarField, plot_kymographs, MemoryStorage, FileStorage
from pde.solvers import Controller, ExplicitMPISolver
from pde.tools import mpi 

#defining the PDE class object for the BWH system containing the paramters
#B = vegetation, W = soil water, H = herbivore
class BWH(PDEBase):             
    def __init__(self, P = 120, bc = "auto_periodic_neumann"):
        super().__init__
        self.lamda = 0.5            #vegetation growth rate
        self.E = 10.0               #root to shoot ratio
        self.K = 0.9                #maximum vegetation per unit area
        self.M = 11.4               #vegetation mortality rate
        self.DB = 1.2

        self.N = 20.0               #evaporation rate
        self.R = 0.01               #reduced evaporation due to vegetation
        self.GamW = 10.0            #soil water uptake rate
        self.DW = 150.0             #lateral soil water diffusion coefficient
        
        self.mH = 0.06              #herbivore mortality rate
        self.GamH = 0.3             #herbivore prolifereation rate
        self.AlphaH = 0.6           #vegetation consumption rate
        self.BetaH = 0.82           #satiation biomass
        self.K_H = 150.0            #maximum herbivore per unit area
                    
        self.DHH = 400              #diff. coeff. for random motility
        self.Zeta = 0.001           #reference biomass for 50% random motility drop
        self.DHB = 700              #diff. coeff. for vegetaxis motility
        self.k = 0.0001             #reference biomass for 50% vegetaxis motitity drop
        
        self.P = P                  #precipitation
        self.bc = bc                #boundary condition
    

    #definition of the function for initial condition
    #if a special initial condition is used it can be input through Bini, Wini, Hini during call
    #def get_initial_state(self, grid, storage, Bini, Wini, Hini):
    
    #if random initial condition is used then no need input of Bini, Wini, Hini
    def get_initial_state(self, grid, storage, Lx, mesh):
        # Random initial condition
        B = ScalarField.random_uniform(grid)
        W = ScalarField.random_uniform(grid)
        H = ScalarField.random_uniform(grid)

        #if initial condition is input during call, modify B, W, H from Bini, Wini, Hini respectively
        #one can also modify the random initial condition by accessing B, W, H data as below
        #B.data = Bini
        #W.data = Wini 
        #H.data = Hini
        
        #example of changing initial condition 
        x = np.linspace(0, Lx, mesh)
        B.data = 0.5 + (0.01 * np.sin(2 * np.pi * 4 * x/Lx))
        W.data = 1.2
        H.data = 0.25

        #this labels are useful when output data file is analysed later
        B.label = "Plants"
        W.label = "Water"
        H.label = "Herbivore"

 
        return FieldCollection([B, W, H])


    #rate equations for B, W, H
    def evolution_rate(self, state, t = 0): 
        B, W, H = state
        

        #rate equations implemented in one dimension
        B_t = (self.lamda * B * W * ((1.0 + (self.E * B))**2) * (1.0 - (B/self.K))) - (self.M * B) - (self.AlphaH * H * B/(self.BetaH + B)) + (self.DB * B.laplace(self.bc))
        W_t = self.P - (self.N * W/(1.0 + (self.R * B/self.K))) - (self.GamW * B * W * (1.0 + (self.E * B))**2) + (self.DW * W.laplace(self.bc))
        H_t = (-self.mH * H) + (self.GamH * self.AlphaH * H * B * (1.0 - (H/self.K_H))/(self.BetaH + B)) - ((self.DHB * self.k/(self.k + B)) * ((H.gradient(self.bc) * B.gradient(self.bc)))) + ((self.DHB * H * self.k/(self.k + B)**2) * B.gradient(self.bc)**2) - ((self.DHB * H * self.k/(self.k + B)) * B.laplace(self.bc)) + ((self.DHH * self.Zeta**2/(self.Zeta**2 + B**2)) * H.laplace(self.bc)) - ((self.DHH * 2 * self.Zeta**2 * B/(self.Zeta**2 + B**2)**2) * ((H.gradient(self.bc) * B.gradient(self.bc))))
        
        #rate equations implemented in two dimension
        #B_t = (self.lamda * B * W * ((1.0 + (self.E * B))**2) * (1.0 - (B/self.K))) - (self.M * B) - (self.AlphaH * H * B/(self.BetaH + B)) + (self.DB * B.laplace(self.bc))
        #W_t = self.P - (self.N * W/(1.0 + (self.R * B/self.K))) - (self.GamW * B * W * (1.0 + (self.E * B))**2) + (self.DW * W.laplace(self.bc))
        #H_t = (-self.mH * H) + (self.GamH * self.AlphaH * H * B * (1.0 - (H/self.K_H))/(self.BetaH + B)) - ((self.DHB * self.k/(self.k + B)) * ((H.gradient(self.bc) @ B.gradient(self.bc)))) + ((self.DHB * H * self.k/(self.k + B)**2) * (B.gradient(self.bc) @ B.gradient(self.bc))) - ((self.DHB * H * self.k/(self.k + B)) * B.laplace(self.bc)) + ((self.DHH * self.Zeta**2/(self.Zeta**2 + B**2)) * H.laplace(self.bc)) - ((self.DHH * 2 * self.Zeta**2 * B/(self.Zeta**2 + B**2)**2) * ((H.gradient(self.bc) @ B.gradient(self.bc))))

        return FieldCollection([B_t, W_t, H_t])
    

    #numba implemention to run this code parallely. Python packages such as Numba, Numba-mpi, mpi4py are required to 
    #run it parallely. We suggest to consult the py-pde documentation page mentioned in the beginning.  
    def _make_pde_rhs_numba(self, state):
        """numba implementation of the PDE"""
        #parameter values as defined before
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
        
        #defining dot operation, required for two dimensional implementation
        #dot = VectorField(state.grid).make_dot_operator()

        #numba implementation of the rate equations
        @nb.jit(nopython=True)
        def pde_rhs(state_data, t):
            B = state_data[0]
            W = state_data[1]
            H = state_data[2]
            
            rate = np.empty_like(state_data)
            
            
            #rate equations implemented in one dimension
            rate[0] = (lamda * B * W * ((1.0 + (E * B))**2) * (1.0 - (B/K))) - (M * B) - (AlphaH * H * B/(BetaH + B)) + (DB * laplace(B))
            rate[1] = P - (N * W/(1.0 + (R * B/K))) - (GamW * B * W * (1.0 + (E * B))**2) + (DW * laplace(W))
            rate[2] = (-mH * H) + (GamH * AlphaH * H * B * (1 - (H/K_H))/(BetaH + B)) - ((DHB * k/(k + B)) * ((gradient(H) * gradient(B)))) + ((DHB * H * k/(k + B)**2) * gradient_squared(B)) - ((DHB * H * k/(k + B)) * laplace(B)) + (((DHH * Zeta**2)/((Zeta**2) + (B**2))) * laplace(H)) - ((DHH * 2 * (Zeta**2) * B/((Zeta**2) + (B**2))**2) * (gradient(H) * gradient(B)))
            
            
            #rate equations implemented in two dimension
            #rate[0] = (lamda * B * W * ((1.0 + (E * B))**2) * (1.0 - (B/K))) - (M * B) - (AlphaH * H * B/(BetaH + B)) + (DB * laplace(B))
            #rate[1] = P - (N * W/(1.0 + (R * B/K))) - (GamW * B * W * (1.0 + (E * B))**2) + (DW * laplace(W))
            #rate[2] = (-mH * H) + (GamH * AlphaH * H * B * (1 - (H/K_H))/(BetaH + B)) - ((DHB * k/(k + B)) * (dot(gradient(H), gradient(B)))) + ((DHB * H * k/(k + B)**2) * gradient_squared(B)) - ((DHB * H * k/(k + B)) * laplace(B)) + (((DHH * Zeta**2)/((Zeta**2) + (B**2))) * laplace(H)) - ((DHH * 2 * (Zeta**2) * B/((Zeta**2) + (B**2))**2) * dot(gradient(H), gradient(B)))
            return rate
        
        return pde_rhs


#domain size in one dimension
#to implement in two dimension, use Ly for domain length in one dimension
Lx = 4*np.pi
mesh = 1024

#defining the grid with periodic boundary condition in one dimension
grid = CartesianGrid([[0,Lx]], [mesh], periodic = True)

#to define the grid with periodic boundary condition in two dimension
#grid = CartesianGrid([[0,Lx],[0,Lx]], [mesh,mesh], periodic = [True,True])

#x = np.linspace(0.0, Lx, mesh)

#___________________________________________________
#reading from data file, the initial condition
#
#f = h5py.File('data.h5','r')   
## This data file is not provided but it can generated in the way we have described in the beginning. 
#keylist = list(f.keys())
#data = f[keylist[0]]
#Bini = data[1000][0]    #here the final time data (frame no: 10000)is read and used as initial condition to Bini and Wini
#Wini = data[1000][1]
#Hini = np.zeros(mesh)   #initialising Hini

#detail of the Hini
#Hini_height = 0.8       #patch height
#Hini_width = 1.2        #patch width
#Hini_position = 22      #patch position, should be decided after plotting Bini and Wini 
#Hini = Hini + (np.exp(-(x - Hini_position)**2/(2. * Hini_width**2)) * Hini_height) #final Hini
#__________________________________________________

P_list = [120]          #precipitation (use 70mm/y for the two dimensional case)
#storage = MemoryStorage()  #useful when using jupyter notebook for the analysis of data without saving a data file   

for p in P_list:
    eq = BWH()          #initialising the PDE class
    eq.P = p            #P input 
    storage = FileStorage("out_data_t30.h5")        #output data file
    trackers = ["progress", storage.tracker(0.5)]   #tracker to see the duration of simulation 

    #defining the solver : Runge-Kutta with adaptive time stepping with domain divided into 4 cores
    solver_mpi = ExplicitMPISolver(eq, scheme = 'rk', decomposition=2, backend = 'auto', adaptive = True)

    #controller for parallel computation specifying the total time duration 
    controller_mpi = Controller(solver_mpi, t_range = 30, tracker = trackers)

    #when Bini, Hini, Wini is provided from data or prepared specially as above then the function call to "get_initial_state" with Bini, Wini, Hini as input 
    #state = eq.get_initial_state(grid, storage, Bini, Wini, Hini)	
    
    #generally if random initial condition is used, then no need for Bini, Wini or Hini as input
    state = eq.get_initial_state(grid, storage, Lx, mesh)
    
    #to run the model !
    sol = controller_mpi.run(state)
