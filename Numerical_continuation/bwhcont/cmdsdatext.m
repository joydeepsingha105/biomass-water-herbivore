%% solution branches
P=zeros(110,3);         %array size depends on the number of points in the dir
for i=1:110;
    pt=['bwh/BS/pt' mat2str(0+i*1)]; tmp=load(pt);tmpp=tmp.p;       %change input directory accordingly
   
    P(i,1)=tmpp.u(tmpp.nu + 6);
    
    np=p.nu/p.nc.neq;n=p.nu;
        
    P(i,2)=(tmpp.u(1:np)'*(p.mat.M(1:np,1:np)*(tmpp.u(1:np))))/tmpp.vol;
    p(i,2)=sqrt(P(i,2));
    
    P(i,3)=tmpp.sol.ineg;
end
ptr=fopen('BS.txt','w');                %change output directory accordingly
fprintf(ptr,'%e\t%e\t%e\n',P');
fclose(ptr);
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