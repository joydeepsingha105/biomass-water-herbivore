function out=hobratw(p,u)  
% hobratw: output for bifurcation diagram, setup for TW-Hopf
% [paras, T,max(u),min(u),||u||]
switch p.sw.spcont
    case 0; np=p.nu/p.nc.neq; n=p.nu; % use nu, i.e. all compos of u 
    case 1; np=p.nu/p.nc.neq/2; n=p.nu/2;
    case 2; np=p.nu/p.nc.neq/2; n=p.nu/2;
    case 3; np=p.nu/p.nc.neq/3; n=p.nu/3;    
end
par=u(p.nu+1:end); % par', pause 
if p.sw.para>2 % Hopf setup; get data from p.hopf
    ho=p.hopf; m1=ho.T; m2=max(max(ho.y(1:n,:))); m3=min(min(ho.y(1:n,:))); 
    l2v=zeros(1,ho.tl); 
    for i=1:length(ho.t); l2v(i)=ho.y(1:n,i)'*(ho.tom.M(1:n,1:n)*ho.y(1:n,i)); end 
    l2=trapz(ho.t,l2v); m4=sqrt(l2/p.vol); 
else % steady state continuation  
    po=getpte(p); try L=p.L; catch; L=po(1,1)-po(1,end); end% domain size 
    m1=L/(p.kwnr*abs(par(p.spar))); % period for TW-cont (1D) 
    m2=max(u(1:n)); 
    %m3=min(u(1:n));
    ppp=getaux(p)';
    m3=ppp(20);
    B_l2=u(1:np)'*(p.mat.M(1:np,1:np)*(u(1:np))); m4=sqrt(B_l2/p.vol);
    W_l2=u(np+1:2*np)'*(p.mat.M(np+1:2*np,np+1:2*np)*(u(np+1:2*np))); m5=sqrt(W_l2/p.vol);
    H_l2=u(2*np+1:3*np)'*(p.mat.M(2*np+1:3*np,2*np+1:3*np)*(u(2*np+1:3*np))); m6=sqrt(H_l2/p.vol);
    
    m7 = sum(u(1:np));
    m8 = sum(u(np+1:2*np));
    m9 = sum(u(2*np+1:3*np)); 
end
out=[par; m1; m2; m3; m4; m5; m6; m7; m8; m9];
    