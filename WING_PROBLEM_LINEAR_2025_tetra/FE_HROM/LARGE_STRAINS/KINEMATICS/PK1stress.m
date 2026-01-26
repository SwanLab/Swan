function PoneST = PK1stress(StwoST,FgradST,ndim) ;
% Computation of PK1 stress tensor from the PK2 stress stacked vector and
% the deformation gradient FgradST
% JAHO, 8-Dec-2020
% See symbolic function /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/FindPK1stress.m
%------------------------------------------
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/
% README_RigidBodyMotions.pdf, page 17
%----------------------------------------------------------------------------------------------------------------------

if ndim == 2
    nF = 4 ;
    nelem_ngaus = size(FgradST,1)/nF  ;
    
    PFCOLS = cell(1,nF) ;
    for icols  =1:nF
        PFCOLS{icols} = icols:nF:size(FgradST,1) ;
    end
    nstrain = 3 ;
    PoneST = zeros(nF*nelem_ngaus,size(FgradST,2))  ;
    SCOLS = cell(1,nstrain) ;
    for icols = 1:nstrain
        SCOLS{icols} = icols:nstrain:size(StwoST,1) ;
    end
    
    %P1 =F1*S1 + F3*S3
    PoneST(PFCOLS{1},:) =  FgradST(PFCOLS{1},:).*StwoST(SCOLS{1},:) + ...
        FgradST(PFCOLS{3},:).*StwoST(SCOLS{3},:)  ;
    %P2 = F2*S2 + F4*S3
    PoneST(PFCOLS{2},:) =  FgradST(PFCOLS{2},:).*StwoST(SCOLS{2},:) + ...
        FgradST(PFCOLS{4},:).*StwoST(SCOLS{3},:)  ;
    
    %P3 = F1*S3 + F3*S2
    PoneST(PFCOLS{3},:) =  FgradST(PFCOLS{1},:).*StwoST(SCOLS{3},:) + ...
        FgradST(PFCOLS{3},:).*StwoST(SCOLS{2},:)  ;
    
    %P4 = F2*S3 + F4*S1
    PoneST(PFCOLS{4},:) =  FgradST(PFCOLS{2},:).*StwoST(SCOLS{3},:) + ...
        FgradST(PFCOLS{4},:).*StwoST(SCOLS{1},:)   ;
    
    
else
    
    nF = 9 ;
    nelem_ngaus = size(FgradST,1)/nF  ;
    
    PFCOLS = cell(1,nF) ;
    for icols  =1:nF
        PFCOLS{icols} = icols:nF:size(FgradST,1) ;
    end
    nstrain = 6 ;
    PoneST = zeros(nF*nelem_ngaus,size(FgradST,2))  ;
    SCOLS = cell(1,nstrain) ;
    for icols = 1:nstrain
        SCOLS{icols} = icols:nstrain:size(StwoST,1) ;
    end
    
    %P1 = F1*S1 + F5*S5 + F6*S6
    PoneST(PFCOLS{1},:) =  FgradST(PFCOLS{1},:).*StwoST(SCOLS{1},:) + ...
        FgradST(PFCOLS{5},:).*StwoST(SCOLS{5},:) + ...
        FgradST(PFCOLS{6},:).*StwoST(SCOLS{6},:)  ;
    %P2 = F2*S2 + F4*S4 + F9*S6
    PoneST(PFCOLS{2},:) =  FgradST(PFCOLS{2},:).*StwoST(SCOLS{2},:) + ...
        FgradST(PFCOLS{4},:).*StwoST(SCOLS{4},:) + ...
        FgradST(PFCOLS{9},:).*StwoST(SCOLS{6},:)  ;
    
    %P3 = F3*S3 + F7*S4 + F8*S5
    PoneST(PFCOLS{3},:) =  FgradST(PFCOLS{3},:).*StwoST(SCOLS{3},:) + ...
        FgradST(PFCOLS{7},:).*StwoST(SCOLS{4},:) + ...
        FgradST(PFCOLS{8},:).*StwoST(SCOLS{5},:)  ;
    
    %P4 = F2*S4 + F4*S3 + F9*S5
    PoneST(PFCOLS{4},:) =  FgradST(PFCOLS{2},:).*StwoST(SCOLS{4},:) + ...
        FgradST(PFCOLS{4},:).*StwoST(SCOLS{3},:) + ...
        FgradST(PFCOLS{9},:).*StwoST(SCOLS{5},:)  ;
    
    %P5 = F1*S5 + F5*S3 + F6*S4
    PoneST(PFCOLS{5},:) =  FgradST(PFCOLS{1},:).*StwoST(SCOLS{5},:) + ...
        FgradST(PFCOLS{5},:).*StwoST(SCOLS{3},:) + ...
        FgradST(PFCOLS{6},:).*StwoST(SCOLS{4},:)  ;
    
    %P6 = F1*S6 + F6*S2 + F5*S4
    PoneST(PFCOLS{6},:) =  FgradST(PFCOLS{1},:).*StwoST(SCOLS{6},:) + ...
        FgradST(PFCOLS{6},:).*StwoST(SCOLS{2},:) + ...
        FgradST(PFCOLS{5},:).*StwoST(SCOLS{4},:)  ;
    
    %P7 = F3*S4 + F7*S2 + F8*S6
    PoneST(PFCOLS{7},:) =  FgradST(PFCOLS{3},:).*StwoST(SCOLS{4},:) + ...
        FgradST(PFCOLS{7},:).*StwoST(SCOLS{2},:) + ...
        FgradST(PFCOLS{8},:).*StwoST(SCOLS{6},:)  ;
    %P8 = F3*S5 + F8*S1 + F7*S6
    PoneST(PFCOLS{8},:) =  FgradST(PFCOLS{3},:).*StwoST(SCOLS{5},:) + ...
        FgradST(PFCOLS{8},:).*StwoST(SCOLS{1},:) + ...
        FgradST(PFCOLS{7},:).*StwoST(SCOLS{6},:)  ;
    %P9 = F2*S6 + F4*S5 + F9*S1
    PoneST(PFCOLS{9},:) =  FgradST(PFCOLS{2},:).*StwoST(SCOLS{6},:) + ...
        FgradST(PFCOLS{4},:).*StwoST(SCOLS{5},:) + ...
        FgradST(PFCOLS{9},:).*StwoST(SCOLS{1},:)  ;
    
    
end