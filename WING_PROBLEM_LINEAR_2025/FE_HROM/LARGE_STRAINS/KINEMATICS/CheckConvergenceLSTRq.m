function [criter_f,norm_res,CONVERGED]=  CheckConvergenceLSTRq(FINT,FEXT,ResF,iter,normd,DATA,qLAT)


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
  if ~isempty(qLAT)
      
  STRING_show = ['iter=%d  resf =',formatSpec,'  resf_ad=',formatSpec,'  res_disp_abs=',formatSpec,' qLAT=',formatSpec,'   \n',] ; 
    fprintf(1,STRING_show,[iter,norm_res,criter_f,normd,qLAT]) ;

  else
      STRING_show =  ['iter=%d  resf =',formatSpec,'  resf_ad=',formatSpec,'  res_disp_abs=',formatSpec,'\n',] ;
          fprintf(1,STRING_show,[iter,norm_res,criter_f,normd]) ;

  end
 
  CONVERGED = 0 ; 
  
  if DATA.NEWTON_RAPHSON.ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS == 0
      % Default option
      if criter_f <= DATA.NEWTON_RAPHSON.TOL_FORCES_REL ||  normd <= DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS
          CONVERGED = 1;
      end
  elseif DATA.NEWTON_RAPHSON.ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS == 1
      if criter_f <= DATA.NEWTON_RAPHSON.TOL_FORCES_REL &&  normd <= DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS
          CONVERGED = 1;
      end
  else
      error('Option not implemented')
  end