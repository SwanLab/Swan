function [criter_f norm_res]=  CheckConvergence(FINT,FEXT,ResF,iter)


norm_fint = norm(FINT);
norm_fext = norm(FEXT);
norm_res = norm(ResF,inf);
%dbstop('8')
if norm_fext>0
    norm_maxf = max([norm_fint,norm_fext,1e-16])  ;
else
    norm_maxf = 1 ;
end
criter_f = norm_res/norm_maxf ;
criter_f = min(criter_f,norm_res) ;

formatSpec = '%10.5e\n';

disp(['iter=',iter,'; resf =',num2str(norm_res,formatSpec),' ; resf_ad=',num2str(criter_f,formatSpec)]) ;
