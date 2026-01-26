function [INFO_RVE,MESHdom,SNAPdisp,DATA,OTHER_output,OPERFE,MATPRO,Fbody,Ftrac,SNAPdisp_AUX,DATAoffline] ...
    = GetDisplacAndMesh_sevDOM(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline)
%--------------------------------------------------------------------------
% FUNCTION: GetDisplacAndMesh_sevDOM
%
% PURPOSE:
%   Extracts and assembles the displacement snapshots and relevant mesh/domain
%   information from multiple training simulations for both the main and auxiliary
%   subdomains, to be used in multiscale reduced-order modeling.
%
%   This function generalizes `GetDisplacAndMesh_1dom` by allowing the treatment
%   of multiple auxiliary domains (subdomains) and multiple training cases.
%   It reconstructs the snapshot matrices (displacement fields) for the main
%   region of interest (ROI) as well as for specified auxiliary domains using
%   previously stored SVD-compressed data.
%
% INPUTS:
%   - NAME_BASE        : Prefix of the parametric study folder (e.g., 'RVE_param_')
%   - DATAcommon       : Common configuration data (time steps, input parameters,
%                        domain/faces, normalization settings, etc.)
%   - NAMEsnap_base    : Base path to snapshot storage folder (e.g., 'SNAPSHOTS/RVE_param_')
%   - DATAoffline      : Offline configuration with optional fields:
%                          • LABEL_DOMAIN             : ID of main domain
%                          • LABELS_FACES             : boundary IDs
%                          • LABEL_AUXDOMAINS         : IDs of auxiliary domains
%                          • AdditionalTests          : alternative snapshot sets
%                          • ProportionSnapshotsToBeTaken : fractions to subsample
%                          • ReturnFullDisplacementMatrix : flag to return full matrices
%
% OUTPUTS:
%   - INFO_RVE         : Info about the main domain (DOF indices, elements, etc.)
%   - MESHdom          : Mesh structure of the main domain
%   - SNAPdisp         : Cell array {#cases} of normalized displacement snapshots 
%                        for the main domain
%   - DATA             : Updated data structure (includes INFO_RVEaux if used)
%   - OTHER_output     : Additional outputs from training stage (Fbody, Ftrac, etc.)
%   - OPERFE           : Element-level FE operators for the domain
%   - MATPRO           : Material model used in the simulations
%   - Fbody            : Cell array of body forces used in each training case
%   - Ftrac            : Cell array of traction vectors used in each training case
%   - SNAPdisp_AUX     : Cell array {#auxDomains × #cases} of normalized
%                        displacement snapshots for auxiliary domains (if defined)
%   - DATAoffline      : Updated offline configuration including domain mappings
%
% FUNCTIONAL DESCRIPTION:
%   - For each training case:
%       • Loads snapshot data (SVD-compressed) from disk
%       • Reconstructs displacement matrices (U*S*V') for main domain
%       • Normalizes snapshot matrices per project (if required)
%       • Optionally subsamples the snapshot matrices
%       • If auxiliary domains are defined:
%           - Maps nodes and elements between reference and auxiliary domains
%           - Applies the same process to generate corresponding snapshots
%
%   - The mappings between the main and auxiliary domains are done using
%     pseudo-centroids and nearest-neighbor matching (via `knnsearch`).
%   - The function ensures that both nodal and Gauss point conversions
%     are stored for further use in ROM projection or basis translation.
%
% USAGE CONTEXT:
%   Typically used during the offline phase of multiscale ROM methods
%   (e.g., EIFEM) for snapshot extraction before mode computation or ECM.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Created: 23-Sep-2023
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

% Modification of GetDisplacAndMesh_1dom.m to cope with several subdomains, JAHO, 23-Sept-2023
% -----------------------------------
if nargin == 0
    load('tmp.mat')
end

Fbody = [] ; Ftrac = [] ;
DATA.INFO_RVEaux = [] ;
DATA.ConversionINDEX_AUXDOM = [] ;

CASES = 1:size(DATAcommon.INPUTS_PARAMETERS,2) ;

DATAoffline = DefaultField(DATAoffline,'AdditionalTests',[] );
DATAoffline = DefaultField(DATAoffline,'ProportionSnapshotsToBeTaken',[] ); % For additional tests (over 1)
DATAcommon = DefaultField(DATAcommon,'NORMALIZATION_SNAPSHOTS_PER_PROJECT',1); % For additional tests (over 1)


DATAoffline = DefaultField(DATAoffline,'ReturnFullDisplacementMatrix',0 );

if ~isempty(DATAoffline.AdditionalTests) ;
    CASES = 1:(size(DATAcommon.INPUTS_PARAMETERS,2) + length(DATAoffline.AdditionalTests)) ;
end



% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
DATAoffline = DefaultField(DATAoffline,'LABEL_AUXDOMAINS',[]) ;
if ~isempty(DATAoffline.LABEL_AUXDOMAINS)
    SNAPdisp_AUX =cell(length(DATAoffline.LABEL_AUXDOMAINS),length(CASES)) ;
    ConversionINDEX_AUXDOM.NodalDOFS = cell(length(DATAoffline.LABEL_AUXDOMAINS),1) ;
    ConversionINDEX_AUXDOM.ELEMENTS = cell(length(DATAoffline.LABEL_AUXDOMAINS),1) ;
else
    SNAPdisp_AUX = [] ;
end


DATAoffline.proportionSNAPSHOTS = ones(size(CASES)) ;
for iproj = 1:length(CASES)
    
    if iproj <= size(DATAcommon.INPUTS_PARAMETERS,2)
        NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    else
        iloc = iproj-size(DATAcommon.INPUTS_PARAMETERS,2) ;
        NAME_FOLDER = DATAoffline.AdditionalTests{iloc} ;
        if ~isempty(DATAoffline.ProportionSnapshotsToBeTaken)
            DATAoffline.proportionSNAPSHOTS(iproj)= DATAoffline.ProportionSnapshotsToBeTaken(iloc) ;
        end
        
    end
    
    
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    
    
    
    
    if iproj == 1
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'MATPRO','MESH','OPERFE','OTHER_output','Fbody','Ftrac') ;
        %----------------------------------------------------
        DATA = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ;
        % --------------------------------------------------------
        DATAoffline = DefaultField(DATAoffline,'LABEL_DOMAIN',1) ;
        IDOM = DATAoffline.LABEL_DOMAIN; % Domain we wish to study
        DATAoffline = DefaultField(DATAoffline,'LABELS_FACES',[ 1,2,3,4]) ;
        IndexBoundaries = DATAoffline.LABELS_FACES ;% Indexes boundaries corresponding to this domain
        [INFO_RVE,MESHdom,OPERFE,MATPRO] = IndexesREFdomain(IDOM,IndexBoundaries,MESH,DATA,OPERFE,MATPRO) ; % Indexes used for extracting information
        
        % Pseudo-centroid reference domain
        CentREF = sum(MESHdom.COOR,1)/size(MESHdom.COOR,1) ;
        COORref = bsxfun(@plus,MESHdom.COOR',-CentREF')' ;
        
        CentroisELEMref = PseudoCentroidsAll_elements(MESHdom.CN,COORref) ;
        
        
        
        % PSEUDO-CENTROIDS OF ALL THE ELEMENTS OF THE REFERENCE DOMAIN
        
        
        DATAoffline  = DefaultField(DATAoffline,'LABEL_AUXDOMAINS',[]) ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Additional domains to be studied %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % We call it "Auxiliar domains"
        IndexBoundaries_AUX = [];
      %  AUX_DOM = [] ;
         ngaus_STRESS = DATA.MESH.ngaus_STRESS ; % Number of Gauss points integration of stresses
         nstrain = DATA.MESH.nstrain; % Number of strain/stress variables
        if ~isempty(DATAoffline.LABEL_AUXDOMAINS)
            INFO_RVEaux = cell(length(DATAoffline.LABEL_AUXDOMAINS),1) ;
            for iaux = 1:length(DATAoffline.LABEL_AUXDOMAINS)
                IDOM  = DATAoffline.LABEL_AUXDOMAINS(iaux) ;
                [ INFO_RVEaux{iaux}, MESHdomAUXloc,~,~] ...
                    = IndexesREFdomain(IDOM,IndexBoundaries_AUX,MESH,DATA,[],[]) ; % Indexes used for extracting information
                
                %% CORRESPONDANCE WITH REFERENCE DOMAIN
                % ---------------------------------------
                % This is where we have to establish the mapping between
                % the reference mesh (MESHdom) and we call the "auxiliar"
                % domains
                % Pseudo-centroid   domain
                CentLOC = sum(MESHdomAUXloc.COOR,1)/size(MESHdomAUXloc.COOR,1) ;
                COORloc = bsxfun(@plus,MESHdomAUXloc.COOR',-CentLOC')' ;
                
                % Pseudo-centroids ALL elements
                CentroisELEMloc = PseudoCentroidsAll_elements(MESHdomAUXloc.CN,COORloc) ;
                
                [IDX,CHECKN] = knnsearch(COORloc,COORref) ;   % So that COORrve{i}(IDX(j),:) = COORref(j,:)
                [IDXelem,CHECKNelem] = knnsearch(CentroisELEMloc,CentroisELEMref) ;   % So that COORrve{i}(IDX(j),:) = COORref(j,:)
                TOL_loc = 1e-4; 
                if norm(CHECKN)/max(abs(COORloc(:))) < TOL_loc
                    % Conversion indexes for nodal DOFs
                    % -----------------------------------------------------------------------
                    ConversionINDEX_AUXDOM.NodalDOFS{iaux} = small2large(IDX,size(COORloc,2)) ;
                    % -----------------------------------------------------------------------
                    ConversionINDEX_AUXDOM.ELEMENTS{iaux} = IDXelem ;
                    ConversionINDEX_AUXDOM.GaussINDEX_scalarSTR{iaux} = small2large(IDXelem,ngaus_STRESS) ;
                    ConversionINDEX_AUXDOM.GaussINDEX_stress{iaux} = small2large(ConversionINDEX_AUXDOM.GaussINDEX_scalarSTR{iaux},nstrain) ;
%                    
% INFO_RVE.GaussINDEX_scalarSTR = small2large(INFO_RVE.BODYELEM_globNUM,ngaus_STRESS);  % Indexes Gauss points for scalar variables (Left-Hand side)
% nstrain = DATA.MESH.nstrain; % Number of strain/stress variables
% INFO_RVE.GaussINDEX_stress = small2large(INFO_RVE.GaussINDEX_scalarSTR ,nstrain);

                    
                    
                    
                else
                    % No mapping is possible
                    error('No mapping found')
                end
                
                
                % CONVERSION ELEMENTS/gAUSS POINTS COMPONENTS/STRESSES
                % We begin by computing the pseudo centroid of each element
                %
                
                
                
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    
    
    
    
    % load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'OTHER_output')
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SELECTED_DOFS = INFO_RVE.DOFS_globNUM ;
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
        DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U(SELECTED_DOFS,:)',SNAP_cluster.DISP.S)' ;
        %    INCLUDE_V = 0;
        %   if DATAoffline.ReturnFullDisplacementMatrix == 1
        DISP_LOC{iloc} = DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
        %  end
        %
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    nsnapshots = ceil(size(SNAPdisp{iproj},2)*min(DATAoffline.proportionSNAPSHOTS(iproj),1)) ;
    SNAPdisp{iproj} = SNAPdisp{iproj}(:,1:nsnapshots) ;
    % % NORMALIZATION
    
    if DATAcommon.NORMALIZATION_SNAPSHOTS_PER_PROJECT  == 1
        
        SNAPdisp{iproj} = SNAPdisp{iproj}/norm(SNAPdisp{iproj},'fro') ; %
        
    end
    
    % AUXILIAR DOMAINS
    % ---------------------------------------------------------------
    if  ~isempty(SNAPdisp_AUX)
        disp('AUXILIAR DOMAINS')
        
        for  iaux = 1:size(SNAPdisp_AUX,1)
            DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
            SELECTED_DOFS = INFO_RVEaux{iaux}.DOFS_globNUM ;
            SELECTED_DOFS = SELECTED_DOFS(ConversionINDEX_AUXDOM.NodalDOFS{iaux}) ;
            for iloc = 1:length(NAME_SNAP_loc)
                Nameloc = NAME_SNAP_loc{iloc} ;
                load(Nameloc,'SNAP_cluster') ;
                % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
                DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U(SELECTED_DOFS,:)',SNAP_cluster.DISP.S)' ;
                %    INCLUDE_V = 0;
                %   if DATAoffline.ReturnFullDisplacementMatrix == 1
                DISP_LOC{iloc} = DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
                %  end
                %
            end
            SNAPdisp_AUX{iaux,iproj} = cell2mat(DISP_LOC);
            SNAPdisp_AUX{iaux,iproj} = SNAPdisp_AUX{iaux,iproj}(:,1:nsnapshots);
            % % NORMALIZATION
                if DATAcommon.NORMALIZATION_SNAPSHOTS_PER_PROJECT  == 1

            SNAPdisp_AUX{iaux,iproj} = SNAPdisp_AUX{iaux,iproj}/norm(SNAPdisp_AUX{iaux,iproj},'fro') ; %
                end
        end
        
        DATA.INFO_RVEaux = INFO_RVEaux;
        DATA.ConversionINDEX_AUXDOM = ConversionINDEX_AUXDOM;
        
    end
    
end


