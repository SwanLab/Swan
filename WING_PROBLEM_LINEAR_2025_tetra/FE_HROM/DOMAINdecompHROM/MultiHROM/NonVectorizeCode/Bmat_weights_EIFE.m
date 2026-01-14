function [Bmat,WEIGHTSinteg] = Bmat_weights_EIFE(DATA,TRANSF_COORD,EIFEoper_all,Vrot)
% B_mat and weights internal forces EIFE method
% JAHO,11-March-2023
if DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1
    indCHOSEN = TRANSF_COORD.IndexParentDomain ;
    WEIGHTSinteg = TRANSF_COORD.detJe*EIFEoper_all(indCHOSEN).INTforces.weights ;  % CECM weights
    Bmat =  (EIFEoper_all(indCHOSEN).INTforces.BmatRED*EIFEoper_all(indCHOSEN).OPER.HdefINV_PsiDEFfT*Vrot)/TRANSF_COORD.SCALEFACTOR;
else
    error('Option not implemented yet')
    % Future implementations (11-March-2023) should consider the
    % preliminary version developed
    %  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/02_POST_PROCESS/NonVectorElastCode/ComputeKeMatrix_DilRot_multi.m
end