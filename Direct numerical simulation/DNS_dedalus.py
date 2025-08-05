# This is the implementation of the biomass-water-herbivore model in Dedalus3. In order to get familer with the dedalus package please visit - https://dedalus-project.org
# This minimal version in dedalus implentation of the model is only used initially to validate the results from using py-pde as this implementation runs slow 
# None of the figures in the manuscript were obtained using this
# For all purposes we reccommend to use DNS_pypde.py and modify that as per need 
import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)

#Parameter

P = 290.0           #precipitation
Lamda = 0.56		#vegetation growth rate
E = 10.0            #root to shoot ratio
K = 0.9             #maximal vegetation per unit square
M = 10.0            #vegetation mortality rate
DB = 1.2            #vegetation diffusion rate as a result of seed dispersal
N = 16.0            #evaporation rate
R = 0.1             #reduced evaporation rate due to vegetation
GamW = 10.0         #water uptake rate
DW = 150.0          #lateral soil water diffusion rate

mH = 0.06           #herbivore mortality rate   
GamH = 0.3          #herbivore proliferation rate
AlphaH = 0.6        #herbivore consumption rate
BetaH = 0.9         #satiation biomass
K_H = 150           #maximal herbivore biomass per unit area
    
DHH = 1000          #diffusion coefficient for random motility 
Zeta = 0.001        #threshold biomass when random motility drops to 50%
    
DHB = 2000          #diffusion coefficient for vegetaxis motility
k = 0.0001          #threshold biomass when vegetaxis motility drops to 50%

Lx = 80             #domain size
Nx = 1024           #grid size
dealias = 3/2       #alising
stop_sim_time = 300 #max time for simulation
timestepper = d3.SBDF2  #time stepper 
timestep = 1e-5         #dt
dtype = np.float64      #variable type


#Bases (Fourier basis on one dimension)
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord,dtype=dtype)
xbasis = d3.RealFourier(xcoord,size=Nx,bounds=(0,Lx),dealias=dealias)

#Fields
b = dist.Field(name='b',bases=xbasis)   #vegetation
w = dist.Field(name='w',bases=xbasis)   #soil water
h = dist.Field(name='h',bases=xbasis)   #herbivore
#Substitutions
dx = lambda A: d3.Differentiate(A, xcoord)  #derivative

#Problem 
problem = d3.IVP([b, w, h], namespace=locals())  #rate equations 
problem.add_equation("dt(b)-(DB*dx(dx(b)))=(Lamda*b*w*((1+(E*b))**2)*(1-(b/K)))-(M*b)-(AlphaH*h*b/(BetaH+b))")                          #vegetation
problem.add_equation("dt(w)-(DW*dx(dx(w)))=P-(N*w/(1+(R*b/K)))-(GamW*b*w*(1+(E*b))**2)")                                                #soil water 
problem.add_equation("dt(h)+(mH*h)= (GamH*AlphaH*h*b/(BetaH+b))-dx((DHB*h*k*dx(b)/(k+b))-(DHH*Zeta*Zeta*dx(h)/((Zeta*Zeta)+(b*b))))")   #herbivores

#Random initial conditions, change for different initial condition  
x = dist.local_grid(xbasis)
b.fill_random('g',seed=42,distribution='uniform')           #vegetation
w.fill_random('g',seed=42,distribution='uniform')           #soil water
h.fill_random('g',seed=42,distribution='uniform')           #herbivore

    
#Solver
solver = problem.build_solver(timestepper)                  #building the solver
solver.stop_sim_time = stop_sim_time                       

count = 0                                                   #counter to keep track of number of steps

while solver.proceed:                                       #deploying the solver 
    solver.step(timestep)
    if solver.iteration % 100000 == 0:
        logger.info('Iteration=%i, Time=%e, dt=%e' %(solver.iteration, solver.sim_time, timestep))
        count = count+1
        #plotting each step 
        l = np.linspace(0, Lx,b['g'].shape[0])
        fig = plt.figure()
        gs = fig.add_gridspec(3, hspace=0)
        axs = gs.subplots(sharex=True, sharey=False)
        fig.suptitle("t = %.2f"%(solver.sim_time))
        axs[0].plot(l, b['g'], color='g', label='plant')
        axs[1].plot(l, w['g'], color='b', label='water') 
        axs[2].plot(l, h['g'], color='r', label='herbivore')
        axs[0].set_ylabel('plant')
        axs[1].set_ylabel('water')
        axs[2].set_ylabel('herbivore')
        axs[2].set_xlabel('x[m]')
        for ax in axs:
            ax.label_outer()
                
        plt.savefig('bwh_%d.png'%(count)        
