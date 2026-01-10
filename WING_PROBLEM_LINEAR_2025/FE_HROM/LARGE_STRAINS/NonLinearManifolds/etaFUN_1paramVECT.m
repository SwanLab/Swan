function [eta,etaDER1,etaDER2] = etaFUN_1paramVECT(qL,DATA)
if nargin == 0
    load('tmp1.mat')
end
%--------------------------------------------------------------------------
% etaFUN_1paramVECT  is an adaptation of tauFUN_1paramVECT,  
[f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL) ;

 
eta = [DATA.UleftSingular*(DATA.SSingular.*f)]' ;
etaDER1 = [DATA.UleftSingular*(DATA.SSingular.*df)]  ;
etaDER2 = [DATA.UleftSingular*(DATA.SSingular.*d2f)]  ;



