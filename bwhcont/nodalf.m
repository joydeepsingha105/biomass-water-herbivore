function f=nodalf(p,u) % BWH model, terms without spatial derivatives 
par=u(p.nu+1:end); up=u(1:p.nu); % split into pars and u 
uf=p.mat.fill*up; % fill to full domain 
n=p.np; B=uf(1:n); W=uf(n+1:2*n); H=uf(2*n+1:3*n); % the 3 fields B, W and H 
Lam=par(1); E=par(2); K=par(3); M=par(4); pp=par(6); N=par(7); R=par(8); 
GamW=par(9); mH=par(11); GamH=par(12); alH=par(13); bH=par(14); KH=par(19);
f1=Lam*B.*W.*(1+E*B).^2.*(1-B/K)-M*B-alH*B.*H./(bH+B); 
f2=pp-N*W./(1+R*B/K)-GamW*B.*W.*(1+E*B).^2;
f3=-mH*H+GamH*alH*B.*H.*(1-(H/KH))./(bH+B);
f=[f1; f2; f3];        % the rhs 
