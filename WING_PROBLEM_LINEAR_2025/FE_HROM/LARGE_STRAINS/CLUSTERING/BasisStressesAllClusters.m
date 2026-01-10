function [BasisStwo_cluster,BasisPone_cluster] = ...
    BasisStressesAllClusters(DATAoffline,SNAPstressSTWO,IND_CLUSTERS,DISP_CONDITIONS,BasisU_cluster,SNAPdisp,OPERFE,MATPRO,DATA,...
    SNAPSHOTS_PER_CLUSTER_after_overlap,BasisS_allclusters)

if nargin == 0
    load('tmp.mat')
end

nclusters = length(SNAPSHOTS_PER_CLUSTER_after_overlap) ;


STRESS_PK2_error = zeros(length(nclusters),1) ;

BasisStwo_cluster = cell(1,length(nclusters)) ;
BasisPone_cluster = cell(1,length(nclusters)) ;

DOFl = DISP_CONDITIONS.DOFl ;
DOFr = DISP_CONDITIONS.DOFr ;
SingularValuesStwo_cluster =  cell(1,length(nclusters)) ;

for icluster = 1:nclusters
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STRESSES DUE TO THE DISPLACEMENTS PROJECTED ON THE LOCAL BASIS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %IND_LOCAL_CLUSTER = find(IND_CLUSTERS == icluster) ;
    
    IND_LOCAL_CLUSTER = SNAPSHOTS_PER_CLUSTER_after_overlap{icluster} ;
    
    [SNAPstressSTWOproj_LOC,SNAPstressPonePROJ_LOC] = ...
        StressesFromDisplacemetsLoc(IND_LOCAL_CLUSTER,icluster,BasisU_cluster,SNAPdisp,DOFr,DOFl,OPERFE,MATPRO,DATA)    ;
    
    SNAPstressSTWO_LOC = SNAPstressSTWO(:,IND_LOCAL_CLUSTER) ;
    
    % Check if it meets the error criterion
    STRESS_PK2_error(icluster) = norm(SNAPstressSTWO_LOC-SNAPstressSTWOproj_LOC,'fro')/norm(SNAPstressSTWO_LOC,'fro')  ;
    disp(['CLUSTER = ',num2str(icluster),'; ERROR stress PK2= ',num2str(STRESS_PK2_error(icluster))]);
    if STRESS_PK2_error(icluster) > DATAoffline.errorSTRESS
        %  dbstop('129')
        error('STRESS ERROR CRITERION NOT MET: TRY WITH A LOWER ERROR TOLERANCE FOR DISPLACEMENTS ')
    end
    
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
    
    % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK2 STRESSES (FOR RECONSTRUCTION PURPOSES)
    % -------------------------------------------------------------------------------------------------
    e0 = DATAoffline.errorPK2stress_basis  ;
    
    mu = 0 ; R = 0 ; DATALOC.RELATIVE_SVD = 1;
    DATALOC.HIDE_OUTPUT = 1 ;
    [BasisStwo_cluster{icluster},SS,VV] =   RSVDT(SNAPstressSTWOproj_LOC,e0,mu,R,DATALOC) ;
    disp(['Number of PK2 stress modes =',num2str(size(BasisStwo_cluster{icluster},2))])
    disp('***********************************************************')
    
    SingularValuesStwo_cluster{icluster} =  SS ;
    
    
    % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK1 STRESSES (FOR applying the ECM )
    % ----------------------------------------------------------------------------------------------
    e0 = DATAoffline.errorPK2stress_basis  ;
    
    mu = 0 ; R = 0 ; DATALOC.RELATIVE_SVD = 1;
    DATALOC.HIDE_OUTPUT = 1 ;
    [BasisPone_cluster{icluster},SS,VV] =   RSVDT(SNAPstressPonePROJ_LOC,e0,mu,R,DATALOC) ;
    disp(['Number of PK1 stress modes =',num2str(size(BasisPone_cluster{icluster},2))])
    disp('***********************************************************')
    % disp('***********************************************************')
    % disp(['Number of PK1 stress modes =',num2str(size(BasisPone,2))])
    % disp('***********************************************************')
    % The basis for the PK1 stresses is weighted with the singular values
    
    BasisPone_cluster{icluster} = bsxfun(@times,BasisPone_cluster{icluster}',SS)' ;
    
    
    
end




if isstruct(SNAPdisp)
    % --------------------------
    % We need  a basis matrix for PK2 stresses (global, for all clusters)
    
    if isempty(BasisS_allclusters)
        for icluster  = 1:nclusters
            BasisStwo_cluster{icluster} = bsxfun(@times,BasisStwo_cluster{icluster}',SingularValuesStwo_cluster{icluster})' ;
        end
        
        
        
        % ---------------------------------------------------
        
        TOL_BLOCK = DATAoffline.TOL_SVD_MATRIX_ALL*ones(length(BasisStwo_cluster),1)' ;
        
        
        DATAsvd=[];
        DATAsvd.HIDE_OUTPUT = 0 ;
        
        disp('Computing basis matrix PK2 stresses (by joining all clusters) ')
        
        [BasisS_allclusters,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(BasisStwo_cluster,TOL_BLOCK,DATAsvd) ;
        disp('***********************************************************')
        disp(['Number of PK2 stress modes (global) =',num2str(size(BasisS_allclusters,2))])
        disp('***********************************************************')
        disp(['Singular Values =',num2str(S')])
        % coeff =cell(size(BasisStwo_cluster)) ;
        %stepsPROJECT = cell(size(BasisStwo_cluster)) ;
        for icluster = 1:length(BasisStwo_cluster)
            coeff = BasisS_allclusters'*BasisStwo_cluster{icluster} ;
            coeff = bsxfun(@times,coeff',1./SingularValuesStwo_cluster{icluster})' ;
            BasisStwo_cluster{icluster} = coeff;
            %  ifin = iini + size(coeff{icluster},2) -1 ;
            %  stepsPROJECT{iproj} = iini:ifin ;
            %  iini = ifin + 1;
        end
        
    else
        
          for icluster = 1:length(BasisStwo_cluster)
             BasisStwo_cluster{icluster}  = BasisS_allclusters'*BasisStwo_cluster{icluster} ;
        %    coeff = bsxfun(@times,coeff',1./SingularValuesStwo_cluster{icluster})' ;
           % BasisStwo_cluster{icluster} = coeff;
            %  ifin = iini + size(coeff{icluster},2) -1 ;
            %  stepsPROJECT{iproj} = iini:ifin ;
            %  iini = ifin + 1;
        end
        
        
    end
    
    
    
    
    
    
    BasisStwo_cluster_str.coeff = BasisStwo_cluster ;
    BasisStwo_cluster_str.BASIS = BasisS_allclusters ;
    
    BasisStwo_cluster = BasisStwo_cluster_str ;
end