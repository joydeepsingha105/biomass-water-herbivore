function r=sG(p,u) % 
f=nodalf(p,u); 
np=p.np; par=u(p.nu+1:end); 
dB=par(5); dW=par(10); dhh=par(15); 
zeta=par(16); zeta2=zeta^2; kap=par(18); dhb=par(17); 
up=u(1:p.nu); uf=p.mat.fill*up;
B=uf(1:np); %W=uf(np+1:2*np); 
H=uf((2*np)+1:3*np);
c1=-dhb*kap*H./(kap+B); c2=dhh*zeta2./(zeta2+B.^2); 
[K31,~,~]=p.pdeo.fem.assema(p.pdeo.grid,c1,0,0);
[K33,~,~]=p.pdeo.fem.assema(p.pdeo.grid,c2,0,0);
K31=filltrafo(p,K31); K33=filltrafo(p,K33);
K=p.mat.K;
K=[dB*K 0*K 0*K; 0*K dW*K 0*K;  K31 0*K K33]; 
s=par(20); F=p.mat.M0*f;
r=K*up-F-s*p.mat.Kx*up; 