function [DATAOUT_cecm,DATA ]= CheckingConsistencyIntegration(DATA,xCECM,wCECM,DATAFUN,DATAapprox)


%DATA.DATAREMOVEPOINTS.PORDER = P;




%%% CHECKING THAT EVERYTHING IS CONSISTENT -- -march 2020
VSinv = [] ;
%eval(NAMEDATA) ;
DATALOC = DATA.DATAREMOVEPOINTS ;
%  [PHIk_y dPHIk_y]= ...
%     EvaluateBasisFunctionAtX(xCECM,VSinv,DATALOC,0,1);
DATALOC.xLIM = DATAFUN.xLIM ;
%DATALOC.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA = DATA.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA ;
DATA = DefaultField(DATA,'PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS',0) ;
DATALOC.PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS = DATA.PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS ;
DATALOC.NAME_FUNCTION_TO_INTEGRATE = DATA.DATAFUN.TYPE ; 

    [PHIk_y,~,DATAOUT_cecm ]= EvaluateBasisFunctionAnalytical(xCECM,VSinv,DATALOC,0,1)  ;
 

EXACT_INTEGRAL = DATAapprox.ExactIntegral ; % Xf'*W ;

APPR_INTEGRAL  = PHIk_y'*wCECM ;

errorMINE = norm(EXACT_INTEGRAL-APPR_INTEGRAL)./norm(EXACT_INTEGRAL)*100 ;
disp('-------------------------------------------------------------------------------------------')
disp(['Error using OUR RULE WITH m = ',num2str(length(wCECM)),' points = ',num2str(errorMINE),' %']) ;
disp('-------------------------------------------------------------------------------------------')