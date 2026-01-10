function [Bmat,WEIGHTSinteg,KcELAS] = Bmat_weights_EIFEbub(DATA,TRANSF_COORD,EIFEoper_all,Vrot)
% Modification of Bmat_weights_EIFE.m to consider bubble modes
% JAHO, 3-oct-2023, UPC, Terrassa
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/04_GeneralTheory.mlx
if nargin == 0
    load('tmp1.mat')
end

if DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1
    indCHOSEN = TRANSF_COORD.IndexParentDomain ;
    WEIGHTSinteg = TRANSF_COORD.detJe*EIFEoper_all(indCHOSEN).INTforces.weights ;  % CECM weights
    
    if ~isempty(EIFEoper_all(indCHOSEN).INFO.DOFsBUB)
        % There are bubble modes
        % --------------------------------------------
        nmodesALL = size(EIFEoper_all(indCHOSEN).INTforces.BmatRED,2); % This is the total number of modes
        ndofsBUB = length(EIFEoper_all(indCHOSEN).INFO.DOFsBUB) ;  % This is the number of bubble modes
        indBOUND = 1:(nmodesALL-ndofsBUB) ;   % These are the indexes of the strain modes 
        indBUB = (nmodesALL-ndofsBUB+1):nmodesALL  ; % These are the indexes of the bubble modes
        % See AUX_PDF/AUX_1.mlx . BmatRED is defined in 
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/EIFE_operatorsBUB.
        % and /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/BmatricesEIFE_v2_CECM.m
        % BmatRED_all =  BstFE*PhiDEF ;
%       % Now we interpolate this matrix
        % BmatRED  =  InterpolationGaussVariablesECM(BmatRED_all,CECM_intforces,DATA.MESH.ngaus_STRESS,DATA.MESH.nstrain) ;
        
        % In the above, PhiDEF = [PhiDEF_strain, GammaBUB]
        
        BmatB = (EIFEoper_all(indCHOSEN).INTforces.BmatRED(:,indBOUND)*EIFEoper_all(indCHOSEN).OPER.HdefINV_PsiDEFfT*Vrot)/TRANSF_COORD.SCALEFACTOR;
        BmatBUB = EIFEoper_all(indCHOSEN).INTforces.BmatRED(:,indBUB)/TRANSF_COORD.SCALEFACTOR ;
        
        % Extended B-matrix
        % --------------------------
        ndim = size(TRANSF_COORD.ROTATION,2) ;
        Bmat = ExtendedBmatrixBUB(BmatB,BmatBUB,ndim,ndofsBUB) ;
        
        % Elastic coarse-scale stiffness 
        KcELAS = ElasticCoarseStiff_1elem(EIFEoper_all,indCHOSEN,TRANSF_COORD,ndofsBUB) ; 

        
    else
        Bmat =  (EIFEoper_all(indCHOSEN).INTforces.BmatRED*EIFEoper_all(indCHOSEN).OPER.HdefINV_PsiDEFfT*Vrot)/TRANSF_COORD.SCALEFACTOR;
          % Elastic coarse-scale stiffness 
          ndofsBUB = 0 ; 
        KcELAS = ElasticCoarseStiff_1elem(EIFEoper_all,indCHOSEN,TRANSF_COORD,ndofsBUB) ; 
    end
else
    error('Option not implemented yet')
    % Future implementations (11-March-2023) should consider the
    % preliminary version developed
    %  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/02_POST_PROCESS/NonVectorElastCode/ComputeKeMatrix_DilRot_multi.m
end