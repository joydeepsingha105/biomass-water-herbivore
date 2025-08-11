%%
%the following section contains the pde2path commands for branch point
%continuation, hopf point continuation, fold point continuation. These
%continuations are extremely sensitive to tolerance whenever H is nonzero.
%In those cases it is suggested to check the results with direct numerical
%simulation. In case there is an issue with the convergence we suggest that
%use direct numerical simulation again to determine the position of the
%branch/hopf/fold points within the desired numerical accuracy. 
%% Branch point continuation (BPC) in par(13)=AlphaH of the Turing BP on UV
p=bpcontini('bwh3/UV','bpt3',13,'bwh3/SP0boundary'); plotsol(p); % init and plot 
p.nc.dlammax=0.01;p.nc.lammax=0.608;p.sw.bifcheck=0; p.nc.lammax=35; % HU
p.sw.foldcheck=0; %switch off fold detection 
p.sol.ds=0.001; p.nc.tol=1e-7; p.nc.del=0.001; p.nc.njthreshsp=1e4; 
p.sw.spjac=0; % numjac for jac of extended system 
p.nc.dsmax=0.01; p.nc.dsmin=0.00001; 
p.file.smod=1; mclf(2); p=cont(p,20);
%% plotting 
f=3; mclf(f); c=6; plotbra('bwh3/SP0boundary',f,c,'cl',p2pc('g1'));  
xlabel('\alpha'); ylabel('P'); grid on; box on; 
%% FPC (in gamma_W) 
p=spcontini('bwh3/SP0','fpt1',9,'bwh3/SP0-fpc'); plotsol(p); pause 
p.nc.dlammax=10; p.nc.lammax=200; hucl; 
p.nc.del=1e-2; p.nc.njthresh=1e-1; p.nc.njthreshsp=1e4; 
p.sol.ds=-0.01; p.plot.bpcmp=p.nc.ilam(2); 
p.nc.tol=5e-5; p.sw.spjac=0;
p.file.smod=1; p.sw.bifcheck=0; p.sw.foldcheck=0; 
p.sw.verb=2; p.nc.dsmax=0.1; p=cont(p,50);
%% HPC
p=hpcontini('bwh3/UH','hpt1',13,'bwh3/UH_hpc'); %hopf point cont
huclean(p); plotsol(p); p.nc.dlammax=0.001; p.nc.lammax=0.61; pause 
p.plot.bpcmp=p.nc.ilam(2);
p.nc.del=1e-1;  %perturbation size for finite differences
p.sw.spjac=0; 
p.nc.njthresh=5e-5; p.nc.njthreshsp=1e7; 
p.sol.ds=-0.001; p.nc.tol=5e-5; p.file.smod=1; 
p.sw.bifcheck=0; p.sw.foldcheck=0; 
p.nc.dsmax=0.1; p.nc.dsmin=0.00001; p=cont(p,20);
%% plotting 
f=3; mclf(f); c=6; plotbra('bwh3/UH_hpcont_1',f,c,'cl',p2pc('g1'),'fp',10);  
xlabel('\alpha'); ylabel('P'); grid on; box on; 
%% Exit from BPC and return to standard continuation w.r.t P  
p=bpcontexit('bwh3/SP0boundary','pt10','dummy');  
p=resetc(p); p.sol.ds=-0.0001;p.nc.dsmax=0.001; 
p.nc.tol=1e-6;
p.plot.bpcmp=24;  %switch for selecting cmp for plot
p.sw.spcalc=1;    %switch for calculating stability eigenvalues
clf(2); p.nc.nq=1; p.nc.ilam=[6,20]; 
p.fuha.qf=@qf1; p.fuha.qfder=@qf1der; p.sw.qjac=1; %switch pc on
p.nc.lammax=600;p=cont(p,50);

