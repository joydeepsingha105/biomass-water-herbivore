%% driver script to compute BD for "Set A" with alpha=0.48; run cell-by-cell 
close all; keep pphome;rng(42);
%% specification of the system, parameter values, domain size, grid size 
setApar;  % collects parameters from set A, except AlphaH, which we set here; 
AlphaH=0.603;    % maximum rate of plant consumption per unit herbivore
% parameters are refered to by their number in p.nc.ilam; note numbers here! 
par=[Lambda,E,K,M,DB,P,N,R,GamW,DW,mH,GamH,AlphaH,BetaH,DHH,Zeta,DHB,k,K_H,s];
%       1            6                 12                         17    
nx=256; lx=pi/2;
p=bwinit(lx,nx,par); p=setfn(p,'bwh3/BS');% init and initial folder-name  
% Bare soil solution branch (BS)
p=cont(p,1000);     %continuation with number of steps
%% branch switching to uniform vegetation solution branch without herbivore (UV)
% arguments: input dir, branch point no, output dir, initial stepsize 
p=swibra('bwh3/BS','bpt1','bwh3/UV',-0.05); pause; % pause to inspect kernel in Fig.6
p.nc.dsmax=1.0; p=cont(p,100000);  
%% Uniform solution branch with vegetation and herbivore (UH) using swibra command
p=swibra('bwh3/UV','bpt3','bwh3/UH',-0.01); p.nc.dsmax=5.0; p=cont(p,10000);
%% Extra code (not strictly needed): Uniform solution branch with vegetation and herbivore (UH) 
% through manual initialisation if swibra does'nt work (happens for too large ds); 
p=loadp('bwh3/UV','bpt3','bwh1/UH'); %load from a pt in the in dir to an out dir
p.u(p.nu+6)=p.u(p.nu+6)+1.0; %manually increasing the lam
n=p.nu/3; p.u(2*n+1:3*n)=2.0;p=resetc(p); %changing H and resetting struc count 
p.sol.ds=1.0; p.nc.dsmax=2; p=cont(p,200);
%% stationary periodic solution branch of vegetation without herbivore (SP0)
p=swibra('bwh3/UV','bpt2','bwh3/SP0',0.01); pause; % pause to see kernel (bif.direction)
p.sw.bifcheck=0; p.nc.dsmax=0.1; p=cont(p,20);pause; % some steps without phase cond.
p.nc.nq=1; p.nc.ilam=[6,20]; %phase condition with speed s as secondary parameter 
p.fuha.qf=@qf1; p.fuha.qfder=@qf1der; %standard phase condition and derivative
p=cont(p,5000);
%% stationary periodic solution branch vegetation and herbivore (SPH)
p=swibra('bwh3/SP0','bpt1','bwh3/SPH',0.1); pause; 
p.nc.tol=1e-8; p.nc.dsmax=1.0; p=cont(p,100);
%% travelling solution branch of both vegetation and herbivore (TW at low P)
% DRIFT bif., hence swibra, not twswibra 
p=swibra('bwh3/SPH','bpt1','bwh3/TW_W',0.01); pause; 
p.nc.tol=1e-5;p.nc.dsmax=1.0; p.sw.bifcheck=2; p=cont(p,150);
%% travelling solution branch from UH branch (TW at high P)
kwnr=1.0; aux.z=[1 -1i]; % aux arg z to get on TW branch 
hp='hpt1'; outb='bwh3/TW_H'; % Hopf point number;output directory 
p=twswibrax('bwh3/UH',hp,20,kwnr,outb,aux); %branch switching from Hopf point to TW 
p.sol.ds=0.1; p.nc.dsmax=0.1; % now setting phase condition (PC);   
p.u0(1:p.nu)=p.tau(1:p.nu); % choose a reference profile u0, and use u0_x for PC 
p.u0=p.u0'; p.u0x=p.mat.Kx*p.u0; plotsolu(p,p.u0x,1,1,1); pause 
p.u(1:p.nu)=p.u(1:p.nu)+0.01*p.tau(1:p.nu); % set predictor by hand 
p.nc.tol=1e-6; p.nc.nq=1; p.nc.ilam=[6;20]; p.fuha.qf=@qf1; p.fuha.qfder=@qjac;  pause 
p.nc.dsmax=0.01; p=cont(p,20); par=getaux(p); par(20)    % to check speed
%% plot branches  
f=6; mclf(f); c=24; 
plotbra('bwh3/BS',f,c,'cl',"#841E0C",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh3/UV',f,c,'cl',"#41E101",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh3/UH',f,c,'cl',"#0066FF",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh3/SP0',f,c,'cl',"#029C7F",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh3/TW_H',f,c,'cl',"#D604D1",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh3/SPH',f,c,'cl','#07F3F9','lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh3/TW_W',f,c,'cl',"#D604D1",'lwun',2,'lwst',5,'ms',0,'fms',0);
%% plot solutions
plotsol('bwh3/SPH','pt100'); 
%% continuation wrt DHB
p=loadp('bwh3/TW_H','pt140','bwh3/speed_DHB');
p.sol.ds=10.0; p.nc.ilam=[17,20]; p=resetc(p); plotsol(p); 
p.nc.lammax=6000; p.nc.tol=1e-4; p.sw.bifcheck=0; p.nc.dsmax=10.0; 
p.plot.bpcmp=23;p=cont(p,3); pause; 
p.nc.nq=1;  p.fuha.qf=@qf1; 
p.fuha.qfder=@qf1der; p.sw.qjac=1; % switch on PC 
p.sw.bifcheck=2;p=cont(p,1000);

    