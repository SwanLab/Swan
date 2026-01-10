function [w,wDER1,wDER2] = wFUN_1paramVECTp(qL,DATA)
if nargin == 0
    load('tmp1.mat')
end
%--------------------------------------------------------------------------
% wFUN_1paramVECT  is an adaptation of tauFUN_1paramVECT,  

% evaluate_spline_with_extrapolationMAWecm is the approach using projection
% onto the active set,
% see % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/10_SAW_ECMlarge.mlx
%[w,wDER1,wDER2 ] = evaluate_spline_with_extrapolationMAWecm(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL,DATA.VOL) ;

% This method consists simply in "freezing" the weights outside the
% training region
[w,wDER1,wDER2 ] = evaluate_spline_with_extrapolationMAWecm(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL,DATA.VOL) ;


% w = DATA.VOL*w ; 
% wDER1 = DATA.VOL*wDER1 ; 
% wDER2 = DATA.VOL*wDER2 ; 

 
% w = [DATA.UleftSingular*(DATA.SSingular.*f)]' ;
% wDER1 = [DATA.UleftSingular*(DATA.SSingular.*df)]  ;
% wDER2 = [DATA.UleftSingular*(DATA.SSingular.*d2f)]  ;



