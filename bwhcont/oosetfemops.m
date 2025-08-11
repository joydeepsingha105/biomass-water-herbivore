function p=oosetfemops(p)
gr=p.pdeo.grid; fem=p.pdeo.fem; 
[K,M,~]=fem.assema(gr,1,1,1); % (scalar) stiffness and mass matrices
M=kron([[1,0,0];[0,1,0];[0,0,1]],M); % system mass matrix 
p.mat.M0=p.mat.fill'*M; % adapting to pBCs, here for nonlinearity
p.mat.M=filltrafo(p, M); p.mat.K=filltrafo(p, K); % adapting to pBCs (general)
Kx=fem.convection(gr,1); % advection (scalar), now transform to system and pBCs 
Kx=kron([[1,0,0];[0,1,0];[0,0,1]],Kx); p.mat.Kx=filltrafo(p,Kx);
