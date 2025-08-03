%% solution branches
npt=3; nski=20; P=zeros(npt,3);         %array size depends on the number of points in the dir
for i=1:npt;
    pt=['bwh/BS/pt' mat2str((i-1)*nski)]; tmp=load(pt);tmpp=tmp.p;       %change input directory accordingly
   
    P(i,1)=tmpp.u(tmpp.nu + 6);
    
    np=p.nu/p.nc.neq;n=p.nu;
        
    P(i,2)=(tmpp.u(1:np)'*(p.mat.M(1:np,1:np)*(tmpp.u(1:np))))/tmpp.vol;
    P(i,2)=sqrt(P(i,2));
    
    P(i,3)=tmpp.sol.ineg;
end
ptr=fopen('BS.txt','w');                %change output directory accordingly
fprintf(ptr,'%e\t%e\t%e\n',P');
fclose(ptr);
%% it's all on the branch! (the trick is to find it!); -see also bradat.m 
% branch=[count; type; ineg; lam; err; ||u_1||_2; userdata]
%          1                                      7 .....
% here userdata=output of hobrax=[par; m1; m2; m3; m4;  ...; m9], 
% and m4=normalized-L2 of u1=B is what we want. Since length(par)=20, the
% right index of m1 is 6+20+4=30
% (In plotbra, the first 6 are not counted, hence there m4 is index c=24). 
p=loadp('bwh/BS'); % loads last point from dir 
bra=p.branch; ineg=bra(3,:); lamv=bra(4,:); l2B=bra(30,:); 
P2=[lamv', l2B', ineg']; % this now contains all data from ALL THE COMPUTED points, 
% no matter how many points are saved (no matter what smod is) 
% Of course, this needs that you know in advance what you want. If you need
% something which is not computed in p.fuha.outfu (happens to me!), then
% you must do it your way...
%% check for SPH: 
npt=8; nski=20; P=zeros(npt,3);         %array size depends on the number of points in the dir
for i=1:npt;
    pt=['bwh/SPH/pt' mat2str((i-1)*nski)]; tmp=load(pt);tmpp=tmp.p;          
    P(i,1)=tmpp.u(tmpp.nu + 6); np=p.nu/p.nc.neq;n=p.nu;        
    P(i,2)=(tmpp.u(1:np)'*(p.mat.M(1:np,1:np)*tmpp.u(1:np)))/tmpp.vol; % (normalized) L2 of B
    P(i,2)=sqrt(P(i,2));
    P(i,3)=tmpp.sol.ineg;
end
ptr=fopen('SPH.txt','w');                %change output directory accordingly
fprintf(ptr,'%e\t%e\t%e\n',P'); fclose(ptr);
%% 
p=loadp('bwh/SPH','pt20'); bra=p.branch; ineg=bra(3,:); lamv=bra(4,:); l2B=bra(30,:); 
P2=[lamv', l2B', ineg']; P2
%% branch, fold, hopf points 
P=zeros(27,3);      %array size depends on the number of points in the dir
for i=1:27;
    pt=['bwh/BS' ...
        '/pt' mat2str(0+i*1)]; tmp=load(pt);tmpp=tmp.p; %change input directory accordingly
    P(i,1)=tmpp.u(tmpp.nu + 6);

    P(i,2) = tmpp.u(tmpp.nu + 13);
    P(i,3)=tmpp.sol.ineg;    

end
ptr=fopen('bwh/UV2bpcont.txt','w');        %change output directory accordingly
fprintf(ptr,'%e\t%e\t%e\n',P');
fclose(ptr);
%% solution profile data
n = p.nu/3;
p=loadp('bwhcont/BS','pt198');
B = p.u(1:n);
W = p.u(n+1:2*n);
H = p.u(2*n+1:3*n);
ptr=fopen('bwhcont/pt198.txt','w');
fprintf(ptr,'%f\t%f\t%f\n',[B W H]'); 
fclose(ptr);