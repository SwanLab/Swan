function [QrotNEW,norm_angle_rads] = UpdateRotationsEIFEM2(OPERFE,Qrot,VAR,DATA,iterROTATIONS,dQrot) ;
% Compute incremental rotations, and update rotation matrices (EIFEM, co-rotm v2)
%Latex notation:: /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTexAPPV.tex 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/04_UPDrotUNC.mlx%
% JAHO, 3-dec-2024, THURSDAY, UPC, TERRASSA. 
% ---------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

% \QrotE{e} \leftarrow \dQrotE{e} \QrotE{e} 

 QrotNEW =  MultiplyMatrixBlocks(dQrot,Qrot) ;  
 if DATA.MESH.ndim == 2
     cos_ang = dQrot(1:2:end,1); 
     ang_rads = acos(cos_ang) ; 
     norm_angle_rads = norm(ang_rads)/length(ang_rads) ; 
     
 elseif DATA.MESH.ndim == 3
     
     error('To be implemented')
 else
     error('Option not implemented')
 end 
 
 
 
formatSpec = '%10.5e'; 

 % disp(['iter=',num2str(iter),'; resf =',num2str(norm_res,formatSpec),' ; resf_ad=',num2str(criter_f,formatSpec),' res_disp_abs = ',num2str(normd,formatSpec)]) ;
  
  STRING_show = ['****** iterROT=%d  norm_angle =',formatSpec,'  strain_energy=',formatSpec ,'   \n'] ; 
  fprintf(1,STRING_show,[iterROTATIONS,norm_angle_rads,VAR.STRAIN_ENERGY]) ;