function [criter_f,norm_res,CONVERGED]=  CheckConvergenceLSTR_corot2(FINT,FEXT,ResF,iter,normd,DATA,ResCOMPall)
% Adaptation of CheckConvergenceLSTR.m  to the co-rotational problem,v2
%  JAHO, 02-DEC-2024, UPC, CAMPUS NORD, BARCELONA 
% Latex notation:: /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTexAPPV.tex 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/04_UPDrotUNC.mlx

norm_ResCOMPall = norm(ResCOMPall) ; 
norm_fint = norm(FINT); % Norm of the internal forces
norm_fext = norm(FEXT); % Norm of the external forces 
norm_res = norm(ResF); % Norm of the residual (absolute)
%dbstop('8')
if norm_fext>0
    norm_maxf = max([norm_fint,norm_fext,1e-16])  ;  % We make it dimensionless 
      % But why with respect to the either the internal forces and external forceS?  
else
    norm_maxf = 1 ;
end
criter_f = norm_res/norm_maxf ;


criter_f = min(criter_f,norm_res) ; % What's the goal of this operation ---? We comment it out

formatSpec = '%10.5e'; 

 % disp(['iter=',num2str(iter),'; resf =',num2str(norm_res,formatSpec),' ; resf_ad=',num2str(criter_f,formatSpec),' res_disp_abs = ',num2str(normd,formatSpec)]) ;
  
  STRING_show = ['iter=%d  comp_res =',formatSpec,' resf =',formatSpec,'  resf_ad=',formatSpec,'  res_disp_abs=',formatSpec,'   \n'] ; 
  fprintf(1,STRING_show,[iter,norm_ResCOMPall,norm_res,criter_f,normd]) ;
 
  CONVERGED = 0 ; 
  
  if DATA.NEWTON_RAPHSON.ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS == 0
      % Default option
      if (criter_f <= DATA.NEWTON_RAPHSON.TOL_FORCES_REL ||  normd <= DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS) ...
              && norm_ResCOMPall<=DATA.NEWTON_RAPHSON.TOL_COMPATIBILITY_RESIDUAL
          CONVERGED = 1;
      end
%   elseif DATA.NEWTON_RAPHSON.ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS == 1
%       error('Option not implemented')
% %       if criter_f <= DATA.NEWTON_RAPHSON.TOL_FORCES_REL &&  normd <= DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS
% %           CONVERGED = 1;
% %       end
  else
      error('Option not implemented')
  end