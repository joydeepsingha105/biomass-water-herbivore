%% sample for data extraction for postprocessing (e.g., pyplot) from p.branch 
% branch=[count; type; ineg; lam; err; ||u_1||_2; userdata]
%          1                                      7 .....
% here userdata=output of hobrax=[par; m1; m2; m3; m4;  ...; m9], 
% and m4=normalized-L2 of u1=B is what we want. Since length(par)=20, the
% right index of m1 is 6+20+4=30
% (In plotbra, the first 6 are not counted, hence there m4 is index c=24). 
p=loadp('bwh3/BS'); % loads last point from dir 
bra=p.branch; ineg=bra(3,:); lamv=bra(4,:); l2B=bra(30,:); 
dat=[lamv', l2B', ineg']; 
ptr=fopen('BS.txt','w');         %change output directory accordingly
fprintf(ptr,'%e\t%e\t%e\n',P');   fclose(ptr);
%% alternative: get data from solutions, i.e., loop through branch, load point, 
% compute desired data
npt=3; smod=20; P=[]; 
for i=1:npt;
    pt=['bwh3/BS/pt' mat2str((i-1)*smod)];tmp=load(pt);tmpp=tmp.p; %change input dir as desired
    P(i,1)=tmpp.u(tmpp.nu+6);
    np=p.nu/p.nc.neq;         
    P(i,2)=(tmpp.u(1:np)'*(p.mat.M(1:np,1:np)*(tmpp.u(1:np))))/tmpp.vol;
    P(i,2)=sqrt(P(i,2));
    P(i,3)=tmpp.sol.ineg;
end
ptr=fopen('BS.txt','w');                %change output directory accordingly
fprintf(ptr,'%e\t%e\t%e\n',P'); fclose(ptr);