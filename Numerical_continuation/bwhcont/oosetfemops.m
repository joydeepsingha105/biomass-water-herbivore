function p=oosetfemops(p)
gr=p.pdeo.grid; fem=p.pdeo.fem; 
[K,M,~]=fem.assema(gr,1,1,1); % FEM/mass matrices
M=kron([[1,0,0];[0,1,0];[0,0,1]],M);
p.mat.M0=p.mat.fill'*M;
p.mat.M=filltrafo(p, M);
p.mat.K=filltrafo(p, K);
Kx=fem.convection(gr,1); Kx=kron([[1,0,0];[0,1,0];[0,0,1]],Kx);
p.mat.Kx=filltrafo(p,Kx);
%Dx=makeDx(p); Dx=[Dx,0*Dx,0*Dx; Dx,0*Dx,0*Dx; Dx,0*Dx,0*Dx];
%Dx=filltrafo(p,Dx); n=p.nu/3; p.mat.Dx=Dx(1:n,1:n); 