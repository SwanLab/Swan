function [criter_f, norm_res]=  CheckConvergenceHOMOG(ResF,iter,normd)


norm_res = norm(ResF);
criter_f = norm_res ; 
%dbstop('8')
% if norm_fext>0
%     norm_maxf = max([norm_fint,norm_fext,1e-16])  ;
% % else
% %     norm_maxf = 1 ;
% % end
% criter_f = norm_res/norm_maxf ;
% criter_f = min(criter_f,norm_res) ; 

formatSpec = '%10.5e\n'; 

  disp(['iter=',num2str(iter),'; resf_abs =',num2str(norm_res,formatSpec), ' res_disp_abs = ',num2str(normd,formatSpec)]) ;
  