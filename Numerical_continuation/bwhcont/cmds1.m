%% Please each section of this code to get the branches of 
close all; keep pphome;
%% specification of the system, parameter values, domain size, grid size 
Lambda=0.6;             % rate of vegetation growth per unit soil water
E=10;                   % root to shoot ratio
K=10;                   % maximal standing biomass per unit area
M=11.0;                 % plant mortality rate 
DB=1.2;                 % seed dispersal rate
P=300;                  % precipitation
N=20.0;                 % evaporation rate
R=0.01;                 % Reduction in evaporation due to shading (dimensionless)
GamW=10.0;              % water uptake rate
DW=150;                 % lateral soil water diffusion  
mH=5.0;                 % herbivore mortality rate 
GamH=0.3;               % fraction of consumed biomass used in herbivore production
AlphaH=30;              % maximum rate of plant consumption per unit herbivore
BetaH=0.3;              % satiaton biomass 
DHH=100000;             % maximal random motility
Zeta=0.001;             % reference biomass at which the motility drops to 50%
DHB=5000;               % maximal vegetaxis motility
k=0.0001;               % reference biomass at which the motility drops to 50%
K_H=150;                % maximum herbivore capacity per unit area
s=0;                    % speed
% parameters are refered to by their number in p.nc.ilam; note numbers here! 
par=[Lambda, E, K, M, DB, P, N, R, GamW, DW, mH, GamH, AlphaH, BetaH, DHH, Zeta, DHB, k, K_H, s];
%       1                 6                       12                              17    
nx=256; lx=pi/2;
p=bwinit(lx,nx,par); p=setfn(p,'bwh/BS');% init and initial folder-name  
% Bare soil solution branch (BS)
p=cont(p,1000);     %continuation with number of steps
%% branch switching to uniform vegetation solution branch without herbivore (UV)
% arguments: input dir, branch point no, output dir, initial stepsize 
p=swibra('bwh/BS','bpt1','bwh/UV',0.05); pause; % pause to inspect kernel in Fig.6
p.nc.dsmax=1.0; p=cont(p,100000);  
%% Uniform solution branch with vegetation and herbivore (UH) using swibra command
p=swibra('bwh/UV','bpt2','bwh/UH',-0.01); p.nc.dsmax=5.0;
p.file.smod=20; p=cont(p,10000);
%% Uniform solution branch with vegetation and herbivore (UH) through manual initialisation if swibra does not work 
p=loadp('bwh1/UV','bpt3','bwh1/UH'); %load from a pt in the in dir to an out dir
p.u(p.nu+6)=p.u(p.nu+6)+1.0; %manually increasing the lam
n=p.nu/3; p.u(2*n+1:3*n)=2.0;p=resetc(p); %changing H and resetting struc count 
p.sw.bifcheck=2;        
p.sol.ds=1.0; p.nc.dsmax=2;
p=cont(p,200);
%% stationary periodic solution branch of vegetation without herbivore (SP0)
p=swibra('bwh/UV','bpt4','bwh/SP0',-0.1); pause;
p.sw.bifcheck=0; p.nc.dsmax=0.1; p=cont(p,20);pause;
p.nc.nq=1; p.nc.ilam=[6,20]; %phase condition on and adding speed as another bifurcation parameter 
p.fuha.qf=@qf1; p.fuha.qfder=@qf1der; %standard phase condition and derivative
p.sw.qjac=1;    %phase condition jacobian
p.sw.bifcheck=2;  p=cont(p,5000);
%% stationary periodic solution branch vegetation and herbivore (SPH)
p=swibra('bwh/SP0','bpt3','bwh/SPH',-0.1); pause; 
p.nc.tol=1e-8; p.nc.dsmax=0.1; p=cont(p,200);
%% travelling solution branch of both vegetation and herbivore (TW at low P)
% DRIFT bif., hence swibra, not twswibra 
p=swibra('bwh/SPH','bpt1','bwh/TW',0.1); pause; 
p.nc.tol=1e-5;p.nc.dsmax=1.0; 
p.sw.bifcheck=2; p=cont(p,100000);
%% travelling solution branch from UH branch (TW at high P)
kwnr=1.0;      %guess for spatial wave number
aux.z=[1 -1i]; %auxilary arg for twswibrax, this comb of guesses of z1 and z2 works well for hp point on UH
hp='hpt1'; outb='bwh/TW'; %hopf point number;output directory 
p=twswibrax('bwh/UH',hp,20,kwnr,outb,aux); %branch switching from hopf point to travelling wave branch
p.sol.ds=-0.1; p.nc.dsmax=1;            
p.u0(1:p.nu)=p.tau(1:p.nu);     %reference profile
p.u0=p.u0'; p.u0x=p.mat.Kx*p.u0; % setting PC 
plotsolu(p,p.u0x,1,3,1); 
p.u(1:p.nu)=p.u(1:p.nu)+0.01*p.tau(1:p.nu); 
p.nc.tol=1e-6; p.nc.nq=1; 
p.nc.ilam=[6;20]; p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac;  
p.sw.bifcheck=2; pause;
p.nc.dsmax=0.01; p=cont(p,10); 
par=getaux(p); par(20)          %to check speed

%% continuation wrt DHB
p=loadp('bwh/TW','pt137','bwh/speed_DHB');
p.sol.ds=10.0;p.nc.ilam=[17,20];
p=resetc(p);p.nc.lammax=6000;pause;
p.nc.tol=1e-4; 
p.sw.bifcheck=0; p.nc.dsmax=10.0; 
p.plot.bpcmp=23;p=cont(p,3); pause; 
p.nc.nq=1;  p.fuha.qf=@qf1; 
p.fuha.qfder=@qf1der; p.sw.qjac=1; % switch on PC 
p.sw.bifcheck=2;p=cont(p,1000);

    