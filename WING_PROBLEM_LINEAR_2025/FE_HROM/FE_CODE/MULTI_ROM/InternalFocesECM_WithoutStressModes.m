function  [BasisF, SingVal_F,MSG] = InternalFocesECM_WithoutStressModes(BASES,DATAIN,MSG,BdomRED,Wdom)

if nargin == 0
    load('tmp2.mat')
end

load(DATAIN.NAMEWS_STRESS_SNAPSHOTS,'stressDOM') ;

intfDOM = cell(size(stressDOM));

% Matrix of internal forces snapshots ---directly from stresses
disp('Computing internal forces snapshots')
for itraj = 1:length(stressDOM)
    disp(['Stress block matrix =',num2str(itraj)]) ; 
    USE_SVD_COMPRESS_STRESS = 1;
    if USE_SVD_COMPRESS_STRESS == 1
        % The question here is: should one compress the stresses at this
        % stage ? In principle, it does not make any sense to include all
        % stresses. Why ? Because strains have  already reduced. However,
        % we do not know which are the stresses produced by the alignment
        % method ---nor we need them, don't we. My suggestion (26-Jan-2020) is to compute
        % the internal forces produced, approximately, by truncated stresses. For that, we use the corresponding
        % tolerance, provided by the user. 
        disp('Compressing stresses...')
        %e0=1e-7 ; %,mu,R,DATA
        e0 =DATAIN.TOLERANCE_SVD_RVES_STRESS{itraj}; 
        mu = [] ;
        R = 0 ;
        DATAsvd.RELATIVE_SVD = 1;
        nmodesINI = size(stressDOM{itraj},2) ; 
        [stressDOM{itraj},S,V,eSVD] = RSVDT(stressDOM{itraj},e0,mu,R,DATAsvd ) ;
        stressDOM{itraj} = bsxfun(@times,stressDOM{itraj}',S)' ;
        disp(['From  nsnap = ',num2str(nmodesINI),' to ',num2str(size( stressDOM{itraj} ,2))]) ; 
    end
    
    
    
    intfDOM{itraj} = BasisFfromStress_solo(stressDOM{itraj},...
        BdomRED, Wdom) ;
    % To avoid memory issues, we may  apply the SVD with tolerance 10^-7, to filter
    % noise.
    USE_SVD_COMPRESS=1;
    if USE_SVD_COMPRESS ==1
           disp('Compressing internal forces...')
        e0=1e-7 ; %,mu,R,DATA
        mu = [] ;
        R = 0 ;
        DATAsvd.RELATIVE_SVD = 1;
          nmodesINI = size(intfDOM{itraj},2) ; 
        [intfDOM{itraj},S,V,eSVD] = RSVDT(intfDOM{itraj},e0,mu,R,DATAsvd ) ;
        intfDOM{itraj} = bsxfun(@times,intfDOM{itraj}',S)' ;
        disp(['From  nsnap = ',num2str(nmodesINI),' to ',num2str(size( intfDOM{itraj} ,2))]) ; 
    end
    
    
    
    
end


% Now we apply the SVD to the resulting internal force matrix
nfigure.ERROR = 2789 ;
LEGENDG = 'Internal Forces' ;
DATAIN.NMODES_SHOW = [] ;
RELATIVE_TOLERANCE= DATAIN.TOLERANCE_SVD_RVES_INTERNALF ;
COLOR = 'b';
DATAINloc.LEGEND_GRAPHS = LEGENDG ; 
[BasisF,SingVal_F,V,h1,h2] = SVD_domTOL(intfDOM,nfigure,LEGENDG,DATAIN.NMODES_SHOW,COLOR...
    ,DATAINloc,RELATIVE_TOLERANCE ) ;




%
%      DATAIN.LEGEND_GRAPHS = ['Internal Forces',DATAIN.TypeMode] ;
% nfigure.RIGHTSV = 3200;
% nfigure.ERROR = 1;
% LEGENDG = ['Internal Forces',DATAIN.TypeMode] ;;
% COLOR = 'k' ;
% DATAIN.PARTITION_COLUMNS = DATAIN.nstepsGAUSS ;
% DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_STRESSES ; %(INDICES_PROJECTS) ;
% DATAIN.PLOT_MODES = 1 ;
% DATAIN.ISGAUSS = 1 ;
%
% DATA_REFMESH = [] ;
%
%
% [U,S,V,h1,h2] = SVD_domTOL(dDOM,nfigure,LEGENDG,DATAIN.NMODES_SHOW,COLOR...
%     ,DATAIN,RELATIVE_TOLERANCE ) ;
%
%
% [BasisUdef,SingVal_disp ]= ModesSVDplotUnified(DATAIN,intfDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,RELATIVE_TOLERANCE) ;
%




%
% BasisS = BASES.STRESSES.U ;
% SingVal_stress =  BASES.STRESSES.S/BASES.STRESSES.S(1);  % 28-May-2019
% if length(SingVal_stress) ~= size(BasisS,2)
%     SingVal_stress = ones(size(BasisS,2),1) ;
% end
%
% DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'AlignStressesWithStrains',0) ;
%
% if DATAIN.CUBATURE.AlignStressesWithStrains == 1   && length(SingVal_stress) > size(BdomRED,2)
%     MSG{end+1} = '---------------------------------------------------------' ;
%     MSG{end+1} = 'Alignment with strains when computing integration points' ;
%     MSG{end+1} = ['Number of stress modes before alignment  = ',num2str(length(SingVal_stress))];
%     [BasisS,SingVal_stress] = AlignmentSubspaces(BasisS,SingVal_stress,BdomRED);
% else
%         MSG{end+1} = ['Number of stress modes  = ',num2str(length(SingVal_stress))];
% end
%
% if DATAIN.CUBATURE.IncludeSingularValues_displacements ==1
%    SingVal_Udef =  DATAROM.SingVal_Udef/DATAROM.SingVal_Udef(1) ;
%    BasisUdef = bsxfun(@times,DATAROM.BasisUdef',SingVal_Udef)' ;
%    BdomRED_singval = Bdom*BasisUdef ;
% else
%     BdomRED_singval = BdomRED;
% end
%
%
% % CgloDOM = CgloDOM*Bdom ;
% % Basis for streses
%
% DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'UseStrainsAsStresses',0) ; % Include Singular Value displacement
%
% if DATAIN.CUBATURE.UseStrainsAsStresses == 1
%     BasisS = BdomRED ;
%      SingVal_stress = ones(size(BasisS,2),1) ;
% end
%
%
% DATAINloc.PLOT_RIGHT_SINGULAR_VECTORS  =0  ;
%
% DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'TOL_LOC_InternalForces',1e-6) ;
% DATAINloc.TOL_LOC_InternalForces = DATAIN.CUBATURE.TOL_LOC_InternalForces ;



