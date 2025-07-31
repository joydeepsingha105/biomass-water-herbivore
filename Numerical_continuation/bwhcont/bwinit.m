function p=bwinit(lx,nx,par,varargin)
p=stanparam;  % initialize fields with defaults, then overwrite some: 
p.nc.neq=3; p.sw.sfem=-1; 
p.fuha.outfu=@hobrax; % mod of llibrary function hobra 
p.pdeo=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; % standard 1D PDE object 
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; % # PDE unknowns 
p.nc.neig=30; % #eigenvalues to compute, and reference point for this 
p.sol.xi=1/p.nu; p.file.smod=1; p.sw.para=2; p.sw.foldcheck=1; 
p.nc.ilam=6;  % use par(6) (P) for continuation 
p.sol.ds=0.1; p.nc.dsmax=10; % initial and maximal stepsize 
p.nc.lammin=0.0001; p.nc.lammax=500; % min and max values for cont parameter
p.nc.dlammax=10; %max step size in lam (bifurcation param) 
p.file.smod=20; % save each smod'th cont.step, increase saves less data points 
p.sw.bifcheck=2; % 
p.nc.tol=1e-6; % tolerance 
p.sw.jac=0;    % switch for jacbian, zero as we use numerical jacobian (no sGjac)
p.nc.mu1=1;    % tolerance for entering BP localization by bisection 
%default initial condition
np=p.np; B=zeros(np,1); W=(par(6)/par(7))*ones(np,1); H=zeros(np,1); 
p.u=[B; W; H; par']; % initial solution (bare soil) and parameters
p=box2per(p,1);      % switch on periodic BCs, this also calls oosetfemops  
p.plot.pcmp=[1 2 3]; p.plot.cl={'k','b','r'}; % plotting settings 
screenlayout(p); p.sw.verb=2; % place windows; choose verbosity of output 