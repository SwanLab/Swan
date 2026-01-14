function [z,w,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod(BasisF,SingVal_F,W,DATA)
% See Implementation.pdf
% Empirical Cubature Method

if nargin == 0
    %DATA = [];
    load('tmp.mat')
    DATA.IncludeSingularValuesF = 0 ;
    DATA.TOL = 0 ;
    DATA.EmpiricalCubatureMethod_usingINNERLOOP = 0;
    % G = J ;
end

DATA = DefaultField(DATA,'EmpiricalCubatureMethod_usingINNERLOOP',0)
ERROR_GLO = [] ;  DATAOUT = [] ;
if  DATA.EmpiricalCubatureMethod_usingINNERLOOP == 0
    
    DATA = DefaultField(DATA,'search_in_complementary',0) ; 
    
    if DATA.search_in_complementary == 0
        [z,w,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(BasisF,SingVal_F,W,DATA)  ;% Before 4th 2019. Complete removal of nega. we.
        
    elseif DATA.search_in_complementary == 1
         [z,w,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_CANDcompl(BasisF,SingVal_F,W,DATA)  ;
         
    elseif DATA.search_in_complementary == 2 
        error('This option hasn not proved to offer any advantage ....')

        [z,w,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_CANDcomplRECURS(BasisF,SingVal_F,W,DATA)  ;
    end


else
    error('Option not maintained')
    % It gives essentially the same results as EmpiricalCubatureMethod_orig
    % It was an attempt at using the actual NNLS algorithm,  cf. rutzmoser2018model.pdf, pag 163
    [z,w]= EmpiricalCubatureMethod_NNLS(BasisF,SingVal_F,W,DATA)  ;% After 4th 2019. Sequential removal of nega. we.
end

end