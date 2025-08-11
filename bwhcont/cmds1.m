%% driver script to compute BD for "Set A" with alpha=0.48; run cell-by-cell 
close all; keep pphome; rng(42);
%% specification of the system, parameter values, domain size, grid size 
setApar;  % collect parameters from set A, except AlphaH, which we set now:  
AlphaH=0.48;   % maximum rate of plant consumption per unit herbivore
% parameters are refered to by their number in p.nc.ilam; note numbers here! 
par=[Lambda,E,K,M,DB,P,N,R,GamW,DW,mH,GamH,AlphaH,BetaH,DHH,Zeta,DHB,k,K_H,s];
%       1            6                 12                         17    
nx=256; lx=pi/2;
p=bwinit(lx,nx,par); p=setfn(p,'bwh1/BS');% init and initial folder-name  
% Bare soil solution branch (BS)
p=cont(p,1000);     %continuation with number of steps
%% branch switching to uniform vegetation solution branch without herbivore (UV)
% arguments: input dir, branch point no, output dir, initial stepsize 
p=swibra('bwh1/BS','bpt1','bwh1/UV',-0.05); p.nc.dsmax=1.0; p=cont(p,100000);  
%% Uniform solution branch with vegetation and herbivore (UH) using swibra command
p=swibra('bwh1/UV','bpt3','bwh1/UH',-0.01); p.nc.dsmax=5.0; p=cont(p,10000);
%% Extra code (not strictly needed): Uniform solution branch with vegetation and herbivore (UH) 
% through manual initialisation if swibra does'nt work (happens for too large ds); 
p=loadp('bwh1/UV','bpt3','bwh1/UH'); %load from a pt in the in dir to an out dir
p.u(p.nu+6)=p.u(p.nu+6)+1.0; %manually increasing the lam
n=p.nu/3; p.u(2*n+1:3*n)=2.0;p=resetc(p); %changing H and resetting struc count 
p.sol.ds=1.0; p.nc.dsmax=2; p.file.smod=20; p=cont(p,200);
%% stationary periodic solution branch of vegetation without herbivore (SP0)
p=swibra('bwh1/UV','bpt2','bwh1/SP0',0.01); pause; % pause to see kernel (bif.direction)
p.sw.bifcheck=0; p.nc.dsmax=0.1; p=cont(p,20);pause; % some steps without phase cond.
p.nc.nq=1; p.nc.ilam=[6,20]; %phase condition with speed s as secondary parameter 
p.fuha.qf=@qf1; p.fuha.qfder=@qf1der; %standard phase condition and derivative
p.file.smod=20; p=cont(p,5000);   
%% plotting: branches, see plotbra-source for meaning of arguments  
f=6; mclf(f); c=24; 
plotbra('bwh1/BS',f,c,'cl',"#841E0C",'lwun',2,'lwst',5,'ms',0,'fms',0); pause 
plotbra('bwh1/UV',f,c,'cl',"#41E101",'lwun',2,'lwst',5,'ms',0,'fms',0); pause
plotbra('bwh1/UH',f,c,'cl',"#0066FF",'lwun',2,'lwst',5,'ms',0,'fms',0); pause
plotbra('bwh1/SP0',f,c,'cl',"#029C7F",'lwun',2,'lwst',5,'ms',0,'fms',0); 
xlabel('P'); ylabel('||B||_2'); grid on; box on
%% plot solutions
plotsol('bwh1/SP0','pt100'); 