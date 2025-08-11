%% driver script to compute BD for "Set A" with alpha=0.57; run cell-by-cell 
close all; keep pphome;rng(42);
%% specification of the system, parameter values, domain size, grid size 
setApar;  % collects parameters from set A, except AlphaH, which we set here; 
AlphaH=0.57;    % maximum rate of plant consumption per unit herbivore
% parameters are refered to by their number in p.nc.ilam; note numbers here! 
par=[Lambda,E,K,M,DB,P,N,R,GamW,DW,mH,GamH,AlphaH,BetaH,DHH,Zeta,DHB,k,K_H,s];
%       1            6                 12                         17    
nx=256; lx=pi/2;
p=bwinit(lx,nx,par); p=setfn(p,'bwh2/BS');% init and initial folder-name  
% Bare soil solution branch (BS)
p=cont(p,1000);     %continuation with number of steps
%% branch switching to uniform vegetation solution branch without herbivore (UV)
% arguments: input dir, branch point no, output dir, initial stepsize 
p=swibra('bwh2/BS','bpt1','bwh2/UV',-0.05); p.nc.dsmax=5.0; p.nc.bisecmax=15; p=cont(p,100000);  
%% Uniform solution branch with vegetation and herbivore (UH) using swibra command
p=swibra('bwh2/UV','bpt3','bwh2/UH',-0.01); p.nc.dsmax=5.0; p=cont(p,10000);
%% stationary periodic solution branch of vegetation without herbivore (SP0)
p=swibra('bwh2/UV','bpt2','bwh2/SP0',0.01); pause; % pause to see kernel (bif.direction)
p.sw.bifcheck=0; p.nc.dsmax=0.5; p.nc.tol=1e-4; p=cont(p,4);pause; % some steps without phase cond.
p.nc.tol=1e-6; p.nc.nq=1; p.nc.ilam=[6,20]; %phase condition with speed s as secondary parameter 
p.fuha.qf=@qf1; p.fuha.qfder=@qf1der; %standard phase condition and derivative
p=cont(p,5000);
%% BP4: unphysical!
p=swibra('bwh2/UV','bpt4','bwh2/negH',-0.01); 
%% plot branches  
f=6; mclf(f); c=24; 
plotbra('bwh2/BS',f,c,'cl',"#841E0C",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh2/UV',f,c,'cl',"#41E101",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh2/UH',f,c,'cl',"#0066FF",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh2/SP0',f,c,'cl',"#029C7F",'lwun',2,'lwst',5,'ms',0,'fms',0);
%% plot solutions
plotsol('bwh2/SP0','pt100'); 