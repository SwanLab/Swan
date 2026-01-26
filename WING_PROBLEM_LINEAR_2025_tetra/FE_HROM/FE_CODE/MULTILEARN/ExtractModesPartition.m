function  NAMEWS = ExtractModesPartition(INPUT_DATAFILE,NMODES,DATAIN,DATAONLINE)
% Modal analysis domain-wise,
%dbstop('4')
if nargin ==0
    %%%%%%%%%% Finite Element Projects
    %     INPUT_DATAFILE  = {'DATA_BEAM3DrepMM'} ;    {'DATA_BEAM3DrepMM2', 'DATA_BEAM3DrepMM'}  ;   'DATAStent1' ;  'DATA_BEAMhexaIND'; 'DATA_BEAMhexa'; 'DATA_BEAMthinw';  'DATA_BEAMsolidEMP30'  ;'DATA_BEAMslice'   ;; 'DATA_BEAMsolid' ;{'DATA_BEAMporousBEND'} ; {'DATA_BEAMporousBEND'}; {'DATA_BEAM_Hexa2D_D6_bend'} ; {'DATA_BEAMporousBEND'};;{'DATA_BEAMsolidENA_otherEND'} ; {'DATA_BEAMsolidENA','DATA_BEAMsolidENA_otherEND'} ; {'DATA_BEAM_Hexa2D_D6_bend'} ;  'DATA_BEAM_Hexa2D_bend' ;  % MACRO-STRUCTURE
    %     NMODES= [7] ;   % Number of modes for truncation
    %     DATAIN = [] ;
    load('tmp1.mat')
    
end
if exist(INPUT_DATAFILE{1})~=2
    addpath('DATA_input') ;
end
if exist('ExtractDisplMatrix')~=2
    addpath('FE_CODE') ;
end
if exist('SVDT')==0
    addpath('SVDlibrary')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT inputs
NMODES_SHOW = [] ;
NAMEWS = ['DATAWS/',INPUT_DATAFILE{:},'_',num2str(NMODES),'_OFFLINE','.mat'] ; % Binary file where information is saved

DATAIN = DefaultField(DATAIN,'NAMEWS',NAMEWS);
NAMEWS = DATAIN.NAMEWS;
DATA.NAMEWS = NAMEWS;
DATA.NMODES_TRUNCATE = NMODES;

DATA.TypeUnitCell = 'HEXAG_2D_SQUARE';    % Type of unit cell
DATAIN = DefaultField(DATAIN,'MakeDimensionlessSnapshotsPerProject',0) ;
DIMENSIONLESS = DATAIN.MakeDimensionlessSnapshotsPerProject ;
DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'CUBATURE',0) ; % Efficient integration scheme
DATAIN = DefaultField(DATAIN,'TOLERANCE_SVD_DISPLACEMENTS',[]) ;  % If not-empty, the number of modes are determined by this tolerance
DATAIN = DefaultField(DATAIN,'SEPARATED_MODES',[]) ;  %
DATAIN.SEPARATED_MODES = DefaultField(DATAIN.SEPARATED_MODES,'ACTIVE',0) ;  %  % Separated modes



%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%
%dbstop('42')
if DATAIN.SEPARATED_MODES.ACTIVE == 1
    % Distinct RVEs
    % Number of distinct RVEs
    nTYPErve = length(DATAIN.SEPARATED_MODES.COLUMNS{1}) ;
else
    nTYPErve = 1 ;
end
COLUMNS_RVE = cell(1,nTYPErve) ;  %
COLUMNS_RVEloc = cell(1,length(INPUT_DATAFILE)) ;  %

DATAIN = DefaultField(DATAIN,'REACTIONmodes_equal_DISPmodes',0) ;
if length(DATAIN.REACTIONmodes_equal_DISPmodes) == 1 &  (DATAIN.REACTIONmodes_equal_DISPmodes ~=0)
%    error('The number of components of DATAIN.REACTIONmodes_equal_DISPmodes = number training projects')
end

% -------------------------------------
% CONSTRUCTING MATRICES OF DISP. SNAPSHOTS
% ------------------------------------
iacumPROJ = 0;
if DATAIN.COMPUTE_MODES_AGAIN == 1
    dRVE = [] ;
    %reactDOM = [] ;
    NODESfaces ={} ;
    NODESfaces_ORIG ={} ;
    nDOMprojet = [] ;
    IDXproject = [] ;
    BasisUrb = {};
    disp('--------------------------')
    disp('Collecting disp. snapshots')
    disp('--------------------------')
    MAXd = [] ;
    
    DATAIN = DefaultField(DATAIN,'DOMAINS_TO_INCLUDE_TRAINING',[]) ;
    
    if ~isempty(DATAIN.DOMAINS_TO_INCLUDE_TRAINING)  & DATAIN.SEPARATED_MODES.ACTIVE==1
        error('Non-compatible options')
    end
    selected_columns = [] ;
    for iproject = 1:length(INPUT_DATAFILE)
        disp(['PROJECT = ',INPUT_DATAFILE{iproject}])
        % [dRVEloc,NODESfacesLOC,reactDOMloc ]= ExtractDisplReactionLOC(INPUT_DATAFILE{iproject},DATA)
        %   dbstop('49')
        [dRVEloc,NODESfacesLOC,IDXglo,BasisUrb,NODESfacesLOC_orig,DATAOUT,...
            MaterialTypeLOC,ndim,DOFS_reference ]= ...
            ExtractDisplMatrix(INPUT_DATAFILE{iproject},DATA,DATAIN) ;
        
        if ~isempty(DATAIN.DOMAINS_TO_INCLUDE_TRAINING)
            selected_columnsLOC = DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject} ;
            COLUMNS_RVEloc{iproject} = selected_columnsLOC ; 
            selected_columns = [selected_columns, selected_columnsLOC + iacumPROJ] ;
            %             dRVEloc = dRVEloc(:,selected_columns) ;
            %             NODESfacesLOC = NODESfacesLOC(:,selected_columns) ;
            %             IDXglo = IDXglo(:,selected_columns) ;
            %             NODESfacesLOC_orig = NODESfacesLOC_orig(:,selected_columns) ;
            
        end
        if isempty(COLUMNS_RVEloc{iproject})
        COLUMNS_RVEloc{iproject} = ones(1,size(dRVEloc,2));
        end
        % Determine the columns of dRVEloc corresponding to each type of
        % RVE
        %  dbstop('79')
        [COLUMNS_RVE,COLUMNS_RVEloc ]= ...
            ColumnsSeparatedRVEs(DATAIN,nTYPErve,size(dRVEloc,2),COLUMNS_RVE,iproject,iacumPROJ,COLUMNS_RVEloc);
        %%%%
        iacumPROJ = iacumPROJ + size(dRVEloc,2) ;
        
        %if DIMENSIONLESS == 0
        MAXdLOC = max(sqrt(sum(dRVEloc.*dRVEloc))) ;
        MAXdLOC = ones(1,size(dRVEloc,2))*MAXdLOC ;
        MAXd = [MAXd MAXdLOC] ;
        
        dRVE = [dRVE  (dRVEloc)];
        %else
        %    NORMA   = norm(dRVEloc,'fro') ;
        %    dRVE = [dRVE  (dRVEloc)/NORMA];
        %end
        % reactDOM  = [reactDOM   (reactDOMloc)] ;
        NODESfaces  = [NODESfaces,NODESfacesLOC] ;
        NODESfaces_ORIG  = [NODESfaces_ORIG,NODESfacesLOC_orig] ;
        
        nDOMprojet(iproject) = length(NODESfacesLOC) ;
        IDXproject = [IDXproject IDXglo] ;
        BasisUrbPROJ{iproject} = BasisUrb ; % Rigid body modes
    end
    % dbstop('102')
    if DIMENSIONLESS == 1
        MMM = (max(MAXd)) ; %
        %  dbstop('72')
        dRVE = bsxfun(@times,dRVE',(MMM./MAXd'))' ;
    end
    % % -----------------------------------------------------
    %% Determination of straining modes (SVD of dRVE)
    %------------------------------------------------------
    % dbstop('111')
    if  ~isempty(DATAIN.TOLERANCE_SVD_DISPLACEMENTS) &&  nTYPErve > length(DATAIN.TOLERANCE_SVD_DISPLACEMENTS)
        TOL_LOC = ones(1,nTYPErve)*DATAIN.TOLERANCE_SVD_DISPLACEMENTS(itype) ;
        
    elseif  ~isempty(DATAIN.TOLERANCE_SVD_DISPLACEMENTS)
        TOL_LOC = DATAIN.TOLERANCE_SVD_DISPLACEMENTS;
    else
        TOL_LOC = 0 ;
    end
    BasisUdef ={} ;
    SingVal_disp = {} ;
    
    
    
    
    
    %DATAIN.SVD_based_ON_boundary_DOFS = 1 ;
    DATAIN = DefaultField(DATAIN,'MINIMIZATON_SVD_WITH_STIFFNESS_MATRIX',0) ;
    if DATAIN.MINIMIZATON_SVD_WITH_STIFFNESS_MATRIX == 1
        DATAIN = DefaultField(DATAIN,'SVD_based_ON_boundary_DOFS',1) ;
    else
        DATAIN = DefaultField(DATAIN,'SVD_based_ON_boundary_DOFS',0) ;
    end
    DATAIN.LEGEND_GRAPHS = 'Displacements' ; 
    for itype = 1:nTYPErve
        if ~isempty(TOL_LOC)
            DATAIN.TOL_LOC = TOL_LOC(itype) ;
        else
            DATAIN.TOL_LOC  = [] ;
        end
        nfigure = 1+itype-1;
        LEGENDG = 'Disp. ' ;
        COLOR = 'r' ;
        % dbstop('135')
        
        %         DATAIN = DefaultField(DATAIN,'MINIMIZATON_SVD_WITH_STIFFNESS_MATRIX',0) ;
        %         dbstop('137')
        %         if  DATAIN.MINIMIZATON_SVD_WITH_STIFFNESS_MATRIX == 1 & itype == 1
        %           K = RecoverStiffnessMatrixUnitaryProject(DATAONLINE) ;
        %           Kchol = chol(K) ;
        %           dRVE_loc = Kchol*dRVE(:,COLUMNS_RVE{itype})  ;
        %         else
        %             dRVE_loc = dRVE(:,COLUMNS_RVE{itype}) ;
        %         end
        %
        if ~isempty(selected_columns)
            SEL_COLUMNS = selected_columns ;
        else
            SEL_COLUMNS = COLUMNS_RVE{itype} ;
        end
        DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_DISPLACEMENTS',[]) ; 
        Ui_Si = [] ; 
        if DATAIN.SVD_based_ON_boundary_DOFS == 0
            DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_DISPLACEMENTS ; 
            [U,S,V,h1,h2,Ui_Si] = SVD_and_error(dRVE(:,SEL_COLUMNS),nfigure,LEGENDG,NMODES_SHOW,COLOR...
                ,DATAIN,COLUMNS_RVEloc ) ;
        else
            [U,S,V,h1,h2] = SVD_and_errorBOUNDARY(dRVE(:,SEL_COLUMNS),nfigure,LEGENDG,NMODES_SHOW,...
                COLOR,DATAIN,NODESfaces{1},ndim,DATAONLINE,DOFS_reference,COLUMNS_RVEloc ) ;
        end
        
        if ~isempty(h1)
            hold on
            title(['Domains type =',num2str(itype)])
            legend([h1 ],{'DISPLACEMENTS'})
            legend([h2 ],{'DISPLACEMENTS'})
        end
        
        if isempty(DATA.NMODES_TRUNCATE) & ~isempty(DATAIN.TOL_LOC)
            nnn  =size(U,2) ;
        else
            nnn = min(DATA.NMODES_TRUNCATE,size(U,2)) ;
        end
        BasisUdef{itype} = U(:,1:nnn) ; % = [2] ;
        if isempty(S)
            SingVal_disp{itype} = 0 ;
        else
            SingVal_disp{itype} = S(1:nnn) ;
        end
        
    end
    save(DATA.NAMEWS,'BasisUdef','BasisUrb','NODESfacesLOC','SingVal_disp','COLUMNS_RVE','COLUMNS_RVEloc',...
        'MaterialTypeLOC','DOFS_reference','Ui_Si');
    DATAIN.TOL_LOC =0 ;
    % -----------------------------------------------------------------
    % -------------------------------------
    % CONSTRUCTING MATRICES OF REACTION  SNAPSHOTS (SELF-EQUILIBRATED), AS WELL
    % AS STRESS SNAPSHOTS
    % ------------------------------------
    reactDOM = [] ;
    stressDOM = [] ;
    NODESfaces ={} ;
    iacum  =1;
    disp('--------------------------')
    disp('Collecting reaction  snapshots')
    disp('--------------------------')
    for iproject = 1:length(INPUT_DATAFILE)
        disp(['PROJECT = ',num2str(iproject),' of ',num2str(length(INPUT_DATAFILE)) ])
        % Matrix of snapshot displacements for project "iproject" (referred to reference mesh)
        INDLOC = iacum:iacum+nDOMprojet(iproject)-1 ;
        dRVEloc = dRVE(:,INDLOC) ;
        if ~isempty(IDXproject)
            IDXglo = IDXproject(INDLOC) ; % For referring everything to the original numbering
        else
            IDXglo = [] ;
        end
        NODESfacesLOC = NODESfaces_ORIG(INDLOC) ;
        iacum  = iacum +nDOMprojet(iproject);
        %%%%%
        
        %%%%%
        [reactDOMloc, stressDOMloc,BdomRED, Wdom]= ...
            ExtractReactMatrix(INPUT_DATAFILE{iproject},DATA,dRVEloc,BasisUdef,...
            IDXglo,BasisUrbPROJ{iproject},NODESfacesLOC,DATAIN.CUBATURE,COLUMNS_RVE,...
            COLUMNS_RVEloc{iproject},DATAIN) ;
        
        
        
        
        reactDOM  = [reactDOM   (reactDOMloc)] ;
        %    dbstop('131')
        stressDOM  = [stressDOM   stressDOMloc] ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % SVD of reaction forces
    
    SingVal_reac = {} ; 
    Ui_Si_reac = [] ; 
    for itype = 1:nTYPErve
        
        nfigure = 1+itype-1;
        LEGENDG = 'Reactions. ' ;
        COLOR = 'b' ;
        %dbstop('148')
        if ~isempty(selected_columns)
            SEL_COLUMNS = selected_columns ;
        else
            SEL_COLUMNS = COLUMNS_RVE{itype} ;
        end
        
        %% Rigid body matrix for reactions ---surfaces
        CRIT =  sum(abs(reactDOM),2)  ;  % Determine which components are almost zero
        TOL = max(CRIT)*1e-6;  % Tolerance
        DOF_rest = find(CRIT >=TOL);  % REstricted nodes
        
        DATAIN = DefaultField(DATAIN,'REACTIONmodes_equal_DISPmodes',0);
        if     any(DATAIN.REACTIONmodes_equal_DISPmodes == 1)
            % We ignore the computed value for BasisRdef, and compute it using BasisUdef and BasisUrb
            % First we determine the degrees of freedom in which reactions are
            % zero
            
            
            
            % WE replace the first snapshots of reactions by displacement
            % modes
            NSNAPS = sum(DATAIN.REACTIONmodes_equal_DISPmodes) ;
            % We assume exact projets are placed first in the list ---one
            % spapshot per project
            BASISu = BasisUdef{itype}(:,1:NSNAPS)  ; %  = Deformation modes
            % Self-equilibrated BASIS
            BASIS = BASISu(DOF_rest,:)-  BasisUrb(DOF_rest,:)*(BasisUrb(DOF_rest,:)\BASISu(DOF_rest,:)  ) ;
            reactDOM(:,1:NSNAPS)  = 0 ;
            
            reactDOM(DOF_rest,1:NSNAPS)  = BASIS ;
            
            
        end
        
        DATAIN.no_ortho_first_basis  = 1;
        DATAIN = DefaultField(DATAIN,'REACTION_MODES_FROM_FE_SIMULATIONS',0) ; % = 
        if DATAIN.REACTION_MODES_FROM_FE_SIMULATIONS == 0
            DATAIN.TOLERANCE_SVD_PROJECT =[];
            DATAIN.TOL_LOC = 0 ;
        else
            % 
            
        end
        
        DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_REACTIONS',[]) ;
        
        DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_REACTIONS ;
%         [U,S,V,h1,h2] = SVD_and_error(dRVE(:,SEL_COLUMNS),nfigure,LEGENDG,NMODES_SHOW,COLOR...
%             ,DATAIN,COLUMNS_RVEloc ) ;
%         
       DATAIN.LEGEND_GRAPHS = 'Reactions'; 
        [U,S,V,h1,h2,Ui_Si_reac] = SVD_and_error(reactDOM(:,SEL_COLUMNS),nfigure,LEGENDG,NMODES_SHOW,COLOR,DATAIN,COLUMNS_RVEloc ) ;
        
        if ~isempty(h1)
            hold on
            legend([h1 ],{'REACTIONS'}) ;
            legend([h2 ],{'REACTIONS'}) ;
        end
        nPROJ = sum(DATAIN.REACTIONmodes_equal_DISPmodes) ; 
        if nPROJ >0 & nPROJ ==  size(BasisUdef{itype},2) ;
            nnnnn = size(BasisUdef{itype},2) ;
        else
            nnnnn = size( U,2) ;
        end
        
        % Interface domain method 
        DATAIN = DefaultField(DATAIN,'DOMAIN_DECOMPOSITION_INTERFACE',0) ;  
        DATAIN = DefaultField(DATAIN,'SURPLUS_NUMBER_DISPLACEMENT_MODES',0) ; 
        
        if DATAIN.DOMAIN_DECOMPOSITION_INTERFACE == 1
             nsurplus = DATAIN.SURPLUS_NUMBER_DISPLACEMENT_MODES ;
             nDEFr = size(BasisUdef{itype},2)-nsurplus ; 
             nDEFr = min(nDEFr,size(U,2)) ; 
             BasisRdef{itype} = U(:,1:nDEFr)  ;       
             
             nnnnn = nDEFr ; 
        else
             BasisRdef{itype} = U(:,1:nnnnn)  ; 
        end
        
       
          if isempty(S)
            SingVal_reac{itype} = 0 ;
        else
            SingVal_reac{itype} = S(1:nnnnn) ;
        end
        
    end
    
    %% Rigid body matrix for reactions ---surfaces
    
    BasisRrb = zeros(size(BasisUrb)) ;
    BasisRrb(DOF_rest,:) =  BasisUrb(DOF_rest,:) ;
    
   
    
    
    save(DATA.NAMEWS,'-append','BasisRdef','DATAOUT','stressDOM','BdomRED','Wdom','BasisRrb','SingVal_reac','Ui_Si_reac');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SVD of stresses
    
    
    
else
    % dbstop('170')
    %  if  DATAIN.CUBATURE.ACTIVE == 1
    %     load(DATA.NAMEWS,'BasisS','SingVal_stress','BdomRED','Wdom','DATAOUT','BasisUdef','BasisRdef')
    % else
    load(DATA.NAMEWS,'DATAOUT','BasisUdef','BasisRdef','stressDOM','BdomRED','Wdom','COLUMNS_RVE',...
        'COLUMNS_RVEloc','MaterialTypeLOC','BasisRrb','Ui_Si','Ui_Si_reac')
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduced-order cubature (Empirical Cubature Method)
% ---------------------------------------------------
%dbstop('167')
if   DATAIN.CUBATURE.ACTIVE == 1
    DATAIN.CUBATURE= DefaultField(DATAIN.CUBATURE,'TRUNCATE_FORCES',0) ;
    if DATAIN.COMPUTE_POINTS_AGAIN  == 1
        for itype = 1:nTYPErve
            
            nfigure = 1+itype-1;
            LEGENDG = 'Stresses ' ;
            COLOR = 'k-' ;
            REACTIONmodes_equal_DISPmodes =      DATAIN.REACTIONmodes_equal_DISPmodes ;
            DATAIN.REACTIONmodes_equal_DISPmodes = 0 ;
            [U,S,V,h3,h4] = SVD_and_error(stressDOM(:,COLUMNS_RVE{itype}),...
                nfigure,LEGENDG,NMODES_SHOW,COLOR,DATAIN,COLUMNS_RVEloc ) ;
            
            DATAIN.REACTIONmodes_equal_DISPmodes =   REACTIONmodes_equal_DISPmodes;
            hold on
            legend([h3 ],{'STRESSES'})
            legend([h4 ],{'STRESSES'})
            BasisS{itype} = U  ; % = [2] ;
            SingVal_stress{itype} = S ;
            
        end
        save(DATA.NAMEWS,'-append','BasisS','SingVal_stress');
        % Matrix of internal force snapshot
        % dbstop('228')
        DATAIN.TOL_LOC = 0 ;
        if DATAIN.CUBATURE.TRUNCATE_FORCES == 1
            DATAIN.TOL_LOC = DATAIN.CUBATURE.TOL ;
            DATAIN.CUBATURE.TOL  = 0 ;
        end
        setPoints  = {} ; WdomRED ={} ;
        for itype = 1:nTYPErve
            [BasisF, SingVal_F]= BasisFfromStress(BasisS{itype},SingVal_stress{itype},...
                BdomRED{itype}, Wdom,DATAIN) ;
            
            % Empirical Cubature Method
            DATA.CUBATURE = DATAIN.CUBATURE  ;
            
            [setPoints{itype},WdomRED{itype}]= EmpiricalCubatureMethod(BasisF,SingVal_F,Wdom,DATA.CUBATURE) ;
        end
        save(DATA.NAMEWS,'-append','setPoints','WdomRED')
    else
        load(DATA.NAMEWS,'setPoints','WdomRED')
    end
end

CNref = DATAOUT.CNref;
COOR = DATAOUT.COOR ;
COORref = COOR(DATAOUT.NODESref,:) ;
save(DATA.NAMEWS,'-append','COORref','CNref');

%%% INTERFACE MODES
%  DATAIN = DefaultField(DATAIN,'DOMAIN_DECOMPOSITION_INTERFACE',0) ;  
%  if DATAIN.DOMAIN_DECOMPOSITION_INTERFACE == 1
%      [BasisINT] = BasisInterfaceQ(BasisUdef,SingVal_disp,DATAIN,COORref,CNref,DATA,BasisRdef,DATAOUT) ;
%  else
%      BasisINT = [] ;
%  end


%%%% PLOTTING MODES



%%%%



DATA.PLOT_MODES = 1;



refMESH = 1;
CNref = DATAOUT.CNref; %CNrve{refMESH} ;
NODESref =  DATAOUT.NODESref ;
COOR = DATAOUT.COOR ;
DOFl = [] ;
DATA.NODES = NODESref' ;
TypeElement  =DATAOUT.TypeElement  ;
posgp =  DATAOUT.posgp ;
NAME_MODES = ['DISPLACEMENT'] ;
save(DATA.NAMEWS,'-append','CNref','NODESref','COOR','TypeElement','posgp');

%dbstop('271')
DATA.MaterialType = MaterialTypeLOC ;
if DATA.PLOT_MODES == 1
    
    for itype = 1:length(BasisUdef)
        MODES= BasisUdef{itype} ;
        disp(['Displacements RVE type',num2str(itype),' NUMBER OF MODES = ',num2str(size(BasisUdef{itype},2))])
        
        GidPostProcessModes(COOR,CNref,TypeElement,MODES,posgp,NAME_MODES,DATA,DOFl);
        
        
        %%%% PLOTTING MODES
        disp(['Reactions RVE type',num2str(itype)])
        
        MODES= BasisRdef{itype};
        NAME_MODES = ['REACTIONS'] ;
        %  dbstop('302')
        GidPostProcessModes(COOR,CNref,TypeElement,MODES,posgp,NAME_MODES,DATA,DOFl);
        
    end
    
end
