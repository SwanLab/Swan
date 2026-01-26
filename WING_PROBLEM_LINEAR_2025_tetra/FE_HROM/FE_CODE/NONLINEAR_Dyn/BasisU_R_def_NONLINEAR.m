function [BASES, DATA_REFMESH,nBASES_BEAM,DATAsnap]= BasisU_R_def_NONLINEAR(DATA_TRAINING,DATAIN)
% Computation Basis matrices for displacements, reactions and stresses
% Copy of: BasisU_R_def_computation.m. JAHO, 26-Sept-2018
% -------------------------------------------------------------------
%BASES.DISPLACEMENTS.U =  Basis matrix
%BASES.DISPLACEMENTS.S =  Singular values
%BASES.REACTIONS.U =  Basis matrix
%BASES.REACTIONS.S =  Singular values
%BASES.STRESSES.U =  Basis matrix
%BASES.STRESSES.S =  Singular values
% Likewise, this function return finite element information about the
% reference domain (first domain of the first project), namely

% DATA_REFMESH =
%
%               Nst:  Matrix of stacked shape functions (for Mass  matrix and body forces)
%          WdiagRHS:  Diagonal matrix of weights (for mass matrix and body forces)
%               Bst:  Strain-displacement matrix (global)
%          WdiagLHS:  Diagonal matrix of weights (for stiffness matrix )
%                 K:  Stiffness matrix
%      NameFileMesh:  Name of the mesh file
%     NODES_faces12:  Nodes faces 1 and 2 (sorted)
%            LENGTH:  Length of the domain (x-direction)
%              COOR:  Coordinate matrix
%                CN:  Connectivity matrix
%          CONNECTb{iface}: iface = 1,2 ...nface. Boundary connectivity matrix
%       TypeElement: 'Hexahedra'
%      TypeElementB: 'Quadrilateral'
%             posgp: Gauss points position
%      MaterialType:  Indexes of each material
%              Nbnd{iface}: Matrix of stacked shape functions (for traction forces, face "iface")
%              Wbnd{iface}:  Corresponding matrix of weights
%     DATA_REFMESH.GeometricMassMatrixInterface = Mst ;
%     DATA_REFMESH.RigidBodyModesInterface = R ;


if nargin == 0
    load('tmp1.mat')
end

disp('--------------------------')
disp('Collecting disp. snapshots')
disp('--------------------------')
%  MAXd = [] ;
DATAIN = DefaultField(DATAIN,'DOMAINS_TO_INCLUDE_TRAINING',[]) ;
dDOM = cell(1,length(DATA_TRAINING)) ;  % Displacement snapshots
reactDOM = cell(1,length(DATA_TRAINING)) ;  % REACTIONS snapshots
stressDOM = cell(1,length(DATA_TRAINING)) ;  % REACTIONS snapshots
nstepsNODES = zeros(1,length(DATA_TRAINING))' ;  % 
nstepsGAUSS = zeros(1,length(DATA_TRAINING))' ;  % 

% Check whether some projects are repeated
ISREPEATED = zeros(length(DATA_TRAINING),length(DATA_TRAINING)) ;
nrepeated = zeros(length(DATA_TRAINING),1) ;
for iproj = 1:length(DATA_TRAINING)-1
    for jproj = (iproj+1):length(DATA_TRAINING)
        if strcmp(DATA_TRAINING{iproj},DATA_TRAINING{jproj})
            ISREPEATED(iproj) = jproj  ;
            nrepeated(iproj) = nrepeated(iproj) +1 ;
        end
    end
end

% Boundary interface displacements
DATAIN = DefaultField(DATAIN,'BOUNDARY_INTERFACE_DISPLACEMENT',zeros(length(DATA_TRAINING),1)) ;
DISPLACEMENTS_INTERFACES_BOUNDARY.MODES = [] ;
DISPLACEMENTS_INTERFACES_BOUNDARY.FACES = [] ;

% ----------------------
% Loop over FE projects
% ----------------------
DATA_REFMESH = [] ;
for iproject = 1:length(DATA_TRAINING)
    disp('------------------------------------')
    [~, bbb] = fileparts(DATA_TRAINING{iproject})
    disp(['PROJECT = ',bbb]) ;
    disp('------------------------------------')
    % Reading displacements vectors of each domain (stored in dDOM)
    if isempty(dDOM{iproject})
        [dDOM{iproject},reactDOM{iproject},stressDOM{iproject},DATAREF,nstepsNODES(iproject),nstepsGAUSS(iproject) ]= ...
            ExtractDisplReacMatrix_NONL(DATA_TRAINING{iproject},DATAIN,iproject) ;
    end
    if nrepeated(iproject) > 0
        for irepeated = 1:nrepeated(iproject)
            jproject = ISREPEATED(iproject,irepeated) ;
            dDOM{jproject} = dDOM{iproject} ;
            reactDOM{jproject} = reactDOM{iproject} ;
            stressDOM{jproject} = stressDOM{iproject} ;
        end
    end
    if iproject==1
        DATA_REFMESH = DATAREF ;
    end
    
    %%% Boundary interface displacements
    if DATAIN.BOUNDARY_INTERFACE_DISPLACEMENT(iproject) ==2
        DISPLACEMENTS_INTERFACES_BOUNDARY.MODES = [DISPLACEMENTS_INTERFACES_BOUNDARY.MODES,dDOM{iproject}(:,end)] ;
        DISPLACEMENTS_INTERFACES_BOUNDARY.FACES = [DISPLACEMENTS_INTERFACES_BOUNDARY.FACES; DATAIN.BOUNDARY_INTERFACE_DISPLACEMENT(iproject)] ;
    elseif  DATAIN.BOUNDARY_INTERFACE_DISPLACEMENT(iproject) ==1
        DISPLACEMENTS_INTERFACES_BOUNDARY.MODES = [DISPLACEMENTS_INTERFACES_BOUNDARY.MODES,dDOM{iproject}(:,1)] ;
        DISPLACEMENTS_INTERFACES_BOUNDARY.FACES = [DISPLACEMENTS_INTERFACES_BOUNDARY.FACES; DATAIN.BOUNDARY_INTERFACE_DISPLACEMENT(iproject)] ;
    end
    
    
%     % DomdDOMains to be included in the SVD
%     if ~isempty(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})
%         error('Revise this part of the code !!!! ')
%         selected_columnsLOC = DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject} ;
%         selected_columnsLOC = intersect(1:size(dDOM{iproject},2),selected_columnsLOC) ;
%         dDOM{iproject} = dDOM{iproject}(:,selected_columnsLOC) ;
%         reactDOM{iproject} = reactDOM{iproject}(:,selected_columnsLOC) ;
%         stressDOM{iproject} = stressDOM{iproject}(:,selected_columnsLOC) ;
%     end
end
% % -----------------------------------------------------
%% Determination of BEAM modes (SVD of dDOM)
%------------------------------------------------------
% dbstop('111')
DATAIN = DefaultField(DATAIN,'IS_BEAM_MODE',ones(size(dDOM))) ;


DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_DISPLACEMENTS',cell(size(dDOM))) ;
DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_REACTIONS',cell(size(dDOM))) ;
DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_STRESSES',cell(size(dDOM))) ;

DATAIN = DefaultField(DATAIN,'COMBINED_METHOD',0) ;

DATAIN.nstepsGAUSS = nstepsGAUSS ; 
DATAIN.nstepsNODES = nstepsNODES ; 


if DATAIN.COMBINED_METHOD == 1
    DATAsnap.Displacements =dDOM ;
    DATAsnap.Reactions =reactDOM ;
else
    DATAsnap = [] ;
end

if all(DATAIN.IS_BEAM_MODE)
    % Only beam modes
    ISBEAM=1:length(dDOM) ;
    [BASES.DISPLACEMENTS,BASES.REACTIONS,BASES.STRESSES] = ...
        DeterminationModes(DATAIN,dDOM,reactDOM,...
        stressDOM,DATA_REFMESH,ISBEAM) ;
    
    nBASES_BEAM =[] ;
    
else
    % We have to differentiate between beam modes and enriching modes
    
    ISBEAM = find(DATAIN.IS_BEAM_MODE==1) ; %
    DATAIN.TypeMode = ['_beam'] ;
    [BASES_BEAM.DISPLACEMENTS,BASES_BEAM.REACTIONS,BASES_BEAM.STRESSES] = ...
        DeterminationModes(DATAIN,dDOM(ISBEAM),reactDOM(ISBEAM),...
        stressDOM(ISBEAM),DATA_REFMESH,ISBEAM) ;
    
    nBASES_BEAM.DISPLACEMENTS =  size(BASES_BEAM.DISPLACEMENTS.U,2) ;
    nBASES_BEAM.REACTIONS =  size(BASES_BEAM.REACTIONS.U,2) ;
    nBASES_BEAM.STRESSES =  size(BASES_BEAM.STRESSES.U,2) ;
    close all
    % Now we apply the SVD on the orthogonal complement of the remaining
    % snapshots
    ISNOTBEAM = find(DATAIN.IS_BEAM_MODE==0) ;
    
    dDOM = (dDOM(ISNOTBEAM)) ;
    reactDOM = (reactDOM(ISNOTBEAM)) ;
    stressDOM = (stressDOM(ISNOTBEAM));
    
    MtimesU =DATA_REFMESH.M*BASES_BEAM.DISPLACEMENTS.U ;
    PRB = BASES_BEAM.DISPLACEMENTS.U'*(MtimesU) ;
    
    ORTHOGONAL_WITH_RESPECT_MASS_MATRIX = 1;
    
    
    
    for iprojj = 1:length(dDOM)
        %%%%%%%%
        if ORTHOGONAL_WITH_RESPECT_MASS_MATRIX == 1
            dDOM{iprojj} = dDOM{iprojj} - BASES_BEAM.DISPLACEMENTS.U*(PRB\(MtimesU'*dDOM{iprojj})) ;
            
            %                         if ~isempty(DISPLACEMENTS_INTERFACES_BOUNDARY.MODES)
            %                               DISPLACEMENTS_INTERFACES_BOUNDARY.MODES = DISPLACEMENTS_INTERFACES_BOUNDARY.MODES ...
            %                                   - BASES_BEAM.DISPLACEMENTS.U*PRB*(PRB\(MtimesU'*DISPLACEMENTS_INTERFACES_BOUNDARY.MODES)) ; %(MtimesU\dDOM{iprojj}) ;
            %                         end
            %
            
            
        else
            dDOM{iprojj} = dDOM{iprojj} - BASES_BEAM.DISPLACEMENTS.U*(BASES_BEAM.DISPLACEMENTS.U\dDOM{iprojj}) ;
        end
        
        %%%%%%%%%%%
        reactDOM{iprojj} = reactDOM{iprojj} - BASES_BEAM.REACTIONS.U*(BASES_BEAM.REACTIONS.U\reactDOM{iprojj}) ;
        stressDOM{iprojj} = stressDOM{iprojj} - BASES_BEAM.STRESSES.U*(BASES_BEAM.STRESSES.U\stressDOM{iprojj}) ;
        
    end
    DATAIN.TypeMode = [' localeff'] ;
    DATAIN.NMODES_TRUNCATE = [] ;
    [BASES_ADD.DISPLACEMENTS,BASES_ADD.REACTIONS,BASES_ADD.STRESSES] = ...
        DeterminationModes(DATAIN,dDOM,reactDOM,...
        stressDOM,DATA_REFMESH,ISNOTBEAM) ;
    
    nmodesR  = size(BASES_ADD.REACTIONS.U,2) ;
    nmodesU = size(BASES_ADD.DISPLACEMENTS.U,2) ;
    %nmodesR = min(nmodesR,nmodesU) ; This was eliminated on Sept-2018. 
    
    BASES.DISPLACEMENTS.U = [BASES_BEAM.DISPLACEMENTS.U,  BASES_ADD.DISPLACEMENTS.U] ;
    BASES.DISPLACEMENTS.S = [BASES_BEAM.DISPLACEMENTS.S;  BASES_ADD.DISPLACEMENTS.S ;] ;
    
    
    BASES.REACTIONS.U = [BASES_BEAM.REACTIONS.U,  BASES_ADD.REACTIONS.U(:,1:nmodesR)] ;
    BASES.REACTIONS.S = [BASES_BEAM.REACTIONS.S;  BASES_ADD.REACTIONS.S(1:nmodesR)] ;
    
    %     BASES.REACTIONS.U = [BASES_BEAM.REACTIONS.U,  BASES_ADD.REACTIONS.U] ;
    %     BASES.REACTIONS.S = [BASES_BEAM.REACTIONS.S;  BASES_ADD.REACTIONS.S] ;
    
    BASES.STRESSES.U = [BASES_BEAM.STRESSES.U,  BASES_ADD.STRESSES.U] ;
    BASES.STRESSES.S = [BASES_BEAM.STRESSES.S;  BASES_ADD.STRESSES.S]  ;
    
    
end

BASES.DISPLACEMENTS_INTERFACES_BOUNDARY = DISPLACEMENTS_INTERFACES_BOUNDARY ;

end
