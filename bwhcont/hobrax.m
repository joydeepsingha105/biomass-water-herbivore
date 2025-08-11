function out=hobrax(p,u)  
% hobra: output for bifurcation diagram, setup for Hopf
% [paras, T,max(u),min(u),||u||]
switch p.sw.spcont
    case 0; np=p.nu/p.nc.neq; n=p.nu; % use nu, i.e. all compos of u 
    case 1; np=p.nu/p.nc.neq; n=p.nu/2;
    case 2; np=p.nu/p.nc.neq; n=p.nu/2;
    case 3; np=p.np-1; n=(p.nu-1-p.nc.nq)/3;    
end
%np,n, pause
par=u(p.nu+1:end); % par', pause 
if p.sw.para>2 % Hopf setup; get data from p.hopf
    ho=p.hopf; m1=ho.T; 
    m2=max(max(ho.y(1:np,:))); m3=min(min(ho.y(1:np,:))); 
    l2v=zeros(1,ho.tl); 
    for i=1:length(ho.t); l2v(i)=ho.y(1:np,i)'*(ho.tom.M(1:np,1:np)*ho.y(1:np,i)); end 
    l2=trapz(ho.t,l2v); m4=sqrt(l2/p.vol); 
    m5=0; m6=0; m7=0; m8=0; m9=0; 
else % steady state continuation  
    m1=0; m2=max(u(1:np)); %m3=min(u(1:np));
    
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
    