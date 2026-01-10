function [Nmat,WEIGHTSinteg] = Nmat_weights_EIFE(DATA,TRANSF_COORD,EIFEoper_all)
if nargin == 0
    load('tmp1.mat')
end
% N_mat and weights body forces EIFE method
% JAHO,11-March-2023/18-May-2023
if DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1
    indCHOSEN = TRANSF_COORD.IndexParentDomain ;
    WEIGHTSinteg = TRANSF_COORD.detJe*EIFEoper_all(indCHOSEN).BodyForces.weights ;  % CECM weights
    
    OLD_METHOD_SIMPLE =  0 ;
    if  OLD_METHOD_SIMPLE == 1
        % Before 18-May-2023
        % This is then "old" method ---strictly speaking, it is not
        % correct, although the issue is masked if only the mass matrix is
        % used
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/08_cellullar/05_Vibration_modes.mlx
        Nmat =  EIFEoper_all(indCHOSEN).BodyForces.Nmat ;
    else
        % METHOD AFTER 18-MAY-2023 (STRICTLY CONSISTENT), 
        % see
        % /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/04_Multi2Dcons/EIFEM/EIFEM_pres.pdf
        % SHORTCUT ---it may not be valid for PLATE ELEMENTS !! CAVEAT !!!!
        % 
        Nmat =  EIFEoper_all(indCHOSEN).BodyForces.Nmat ;
        % INTRODUCE ROTATIONS (might not be very efficient, check it)
        % LEFT 
        nDOFSloc =  (size(Nmat,1)) ; 
        nSD = length(TRANSF_COORD.TRANSLATION); 
        nnodesLOC = nDOFSloc/nSD ; 
        ROTglo = cell(1,nnodesLOC) ; 
        ROTglo(:) = {sparse(TRANSF_COORD.ROTATION)} ; 
      %  ROTglo = blkdiag(ROTglo{:}) ; 
        Nmat = blkdiag(ROTglo{:})*Nmat; 
        % RIGHT 
         nDOFSloc =  (size(Nmat,2)) ; 
        nnodesLOC = nDOFSloc/nSD ; 
        ROTglo = cell(1,nnodesLOC) ; 
        ROTglo(:) = {sparse(TRANSF_COORD.ROTATION')} ; 
        Nmat = Nmat*blkdiag(ROTglo{:}); 
    end
else
    error('Option not implemented yet')
    % Future implementations (11-March-2023) should consider the
    % preliminary version developed
    %  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/02_POST_PROCESS/NonVectorElastCode/ComputeKeMatrix_DilRot_multi.m
end