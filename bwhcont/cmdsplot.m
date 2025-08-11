
%% to plot all the solution branchs (mod of plotbra) 
f=6; mclf(f); c=24; 
grid on;
plotbra('bwh1/BS',f,c,'cl',"#841E0C",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh1/UV',f,c,'cl',"#41E101",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh1/UH',f,c,'cl',"#0066FF",'lwun',2,'lwst',5,'ms',0,'fms',0);
plotbra('bwh1/SP0',f,c,'cl',"#029C7F",'lwun',2,'lwst',5,'ms',0,'fms',0);
%plotbra('bwh1/TW_W',f,c,'cl',"#D604D1",'lwun',2,'lwst',5,'ms',0,'fms',0);
%plotbra('bwh1/SPH',f,c,'cl','#07F3F9','lwun',2,'lwst',5,'ms',0,'fms',0);
%plotbra('bwh2/TW_H',f,c,'cl',"#D604D1",'lwun',2,'lwst',5,'ms',0,'fms',0);
xlabel("x[m]");
ylabel("||B||_2[kg/m^2]");  
%% plot solutions
plotsol('bwh/SPH','pt100'); 



