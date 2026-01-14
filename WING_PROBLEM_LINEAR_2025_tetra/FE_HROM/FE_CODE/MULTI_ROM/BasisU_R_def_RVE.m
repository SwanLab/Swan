function [BASES, DATA_REFMESH,nBASES_RVE,DATAsnap,stressDOM,MSG]= BasisU_R_def_RVE(DATA_TRAINING,DATAIN)
% Computation Basis matrices for displacements, reactions and stresses
% ---RVES repeated in the x and y directions
% JAHO, 23 July-2018, copy of BasisU_R_def_computation (for beams)
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

% Disabled, 11-jun-2019
% % Check whether some projects are repeated
% ISREPEATED = zeros(length(DATA_TRAINING),length(DATA_TRAINING)) ;
% nrepeated = zeros(length(DATA_TRAINING),1) ;
% for iproj = 1:length(DATA_TRAINING)-1
%     for jproj = (iproj+1):length(DATA_TRAINING)
%         if strcmp(DATA_TRAINING{iproj},DATA_TRAINING{jproj})
%             ISREPEATED(iproj) = jproj  ;
%             nrepeated(iproj) = nrepeated(iproj) +1 ;
%         end
%     end
% end

% Disabled 11-Jun-2019
% Boundary interface displacements
%DATAIN = DefaultField(DATAIN,'BOUNDARY_INTERFACE_DISPLACEMENT',zeros(length(DATA_TRAINING),1)) ;
%DISPLACEMENTS_INTERFACES_BOUNDARY.MODES = [] ;
%DISPLACEMENTS_INTERFACES_BOUNDARY.FACES = [] ;

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
    %    if isempty(dDOM{iproject})
    [dDOM{iproject},reactDOM{iproject},stressDOM{iproject},...
        DATAREF,nDOMx,nDOMy,timestepINCLUDE,REMOVE_FIRST_DOMAIN  ]= ...
        ExtractDisplReacMatrix_RVE(DATA_TRAINING{iproject},DATAIN,iproject) ;
    %   end
    % disabled, 11-Jun-2019
    %     if nrepeated(iproject) > 0
    %         for irepeated = 1:nrepeated(iproject)
    %             jproject = ISREPEATED(iproject,irepeated) ;
    %             dDOM{jproject} = dDOM{iproject} ;
    %             reactDOM{jproject} = reactDOM{iproject} ;
    %             stressDOM{jproject} = stressDOM{iproject} ;
    %         end
    %     end
    if iproject==1
        DATA_REFMESH = DATAREF ;
    end
    
    %%% Boundary interface displacements
    %     if DATAIN.BOUNDARY_INTERFACE_DISPLACEMENT(iproject) ==2
    %         DISPLACEMENTS_INTERFACES_BOUNDARY.MODES = [DISPLACEMENTS_INTERFACES_BOUNDARY.MODES,dDOM{iproject}(:,end)] ;
    %         DISPLACEMENTS_INTERFACES_BOUNDARY.FACES = [DISPLACEMENTS_INTERFACES_BOUNDARY.FACES; DATAIN.BOUNDARY_INTERFACE_DISPLACEMENT(iproject)] ;
    %     elseif  DATAIN.BOUNDARY_INTERFACE_DISPLACEMENT(iproject) ==1
    %         DISPLACEMENTS_INTERFACES_BOUNDARY.MODES = [DISPLACEMENTS_INTERFACES_BOUNDARY.MODES,dDOM{iproject}(:,1)] ;
    %         DISPLACEMENTS_INTERFACES_BOUNDARY.FACES = [DISPLACEMENTS_INTERFACES_BOUNDARY.FACES; DATAIN.BOUNDARY_INTERFACE_DISPLACEMENT(iproject)] ;
    %     end
    
    
    % DomdDOMains to be included in the SVD
    % Disabled 11-June-2019
    %     if ~isempty(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})
    %
    %         selected_columnsLOC = SelectedComlumnsFun(DATAIN,iproject,nDOMx,nDOMy) ;
    %
    %
    %         selected_columnsLOC = intersect(1:size(dDOM{iproject},2),selected_columnsLOC) ;
    %         dDOM{iproject} = dDOM{iproject}(:,selected_columnsLOC) ;
    %         reactDOM{iproject} = reactDOM{iproject}(:,selected_columnsLOC) ;
    %         stressDOM{iproject} = stressDOM{iproject}(:,selected_columnsLOC) ;
    %     end
    
    
    if REMOVE_FIRST_DOMAIN == 1
        nsteps = length(timestepINCLUDE) ;
        islice = 1;
        COLUMNS_SELECTED = small2large(islice,nsteps) ;
        dDOM{iproject}(:,COLUMNS_SELECTED) = [] ;
        reactDOM{iproject}(:,COLUMNS_SELECTED) = [] ;
        stressDOM{iproject}(:,COLUMNS_SELECTED) = [] ;
    end
    
    
    
end
% % -----------------------------------------------------
%% Determination of RVE modes (SVD of dDOM)
%------------------------------------------------------
% dbstop('111')
DATAIN = DefaultField(DATAIN,'IS_RVE_MODE',ones(size(dDOM))) ;


DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_DISPLACEMENTS',cell(size(dDOM))) ;
DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_REACTIONS',cell(size(dDOM))) ;
DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_STRESSES',cell(size(dDOM))) ;

% DATAIN = DefaultField(DATAIN,'COMBINED_METHOD',0) ;
% if DATAIN.COMBINED_METHOD == 1
%     DATAsnap.Displacements =dDOM ;
%     DATAsnap.Reactions =reactDOM ;
% else
DATAsnap = [] ;
%end

%if all(DATAIN.IS_RVE_MODE)
MSG = {}; 
DATAIN =DefaultField(DATAIN,'CONTROL_TOLERANCE_SVD_PROJECTWISE',0) ; % = 1; 

if DATAIN.ISNONLINEAR == 0  && DATAIN.CONTROL_TOLERANCE_SVD_PROJECTWISE == 0
    % Old method, before 11-June-2019
    ISRVE=1:length(dDOM) ;
    [BASES.DISPLACEMENTS,BASES.REACTIONS,BASES.STRESSES] = ...
        DeterminationModesRVE(DATAIN,dDOM,reactDOM,...
        stressDOM,DATA_REFMESH,ISRVE) ;
    
else
    % Method adapted to nonlinear problems,     %  11-June-2019
    % ------------------------------------

    [BASES.DISPLACEMENTS,BASES.REACTIONS,BASES.STRESSES,MSG] = ...
        DeterminationModesUnified(DATAIN,dDOM,reactDOM,...
        stressDOM,DATA_REFMESH,MSG) ;
end




nBASES_RVE =[] ;

% else  % Option disable 11-June-2019
%     % We have to differentiate between beam modes and enriching modes
%     error('Option not implemented yet')
%     ISRVE = find(DATAIN.IS_RVE_MODE==1) ;
%     DATAIN.TypeMode = ['_beam'] ;
%     [BASES_RVE.DISPLACEMENTS,BASES_RVE.REACTIONS,BASES_RVE.STRESSES] = ...
%         DeterminationModes(DATAIN,dDOM(ISRVE),reactDOM(ISRVE),...
%         stressDOM(ISRVE),DATA_REFMESH,ISRVE) ;
%
%     nBASES_RVE.DISPLACEMENTS =  size(BASES_RVE.DISPLACEMENTS.U,2) ;
%     nBASES_RVE.REACTIONS =  size(BASES_RVE.REACTIONS.U,2) ;
%     nBASES_RVE.STRESSES =  size(BASES_RVE.STRESSES.U,2) ;
%     close all
%     % Now we apply the SVD on the orthogonal complement of the remaining
%     % snapshots
%     ISNOTRVE = find(DATAIN.IS_RVE_MODE==0) ;
%
%     dDOM = (dDOM(ISNOTRVE)) ;
%     reactDOM = (reactDOM(ISNOTRVE)) ;
%     stressDOM = (stressDOM(ISNOTRVE));
%
%     MtimesU =DATA_REFMESH.M*BASES_RVE.DISPLACEMENTS.U ;
%     PRB = BASES_RVE.DISPLACEMENTS.U'*(MtimesU) ;
%
%     ORTHOGONAL_WITH_RESPECT_MASS_MATRIX = 1;
%
%
%
%     for iprojj = 1:length(dDOM)
%         %%%%%%%%
%         if ORTHOGONAL_WITH_RESPECT_MASS_MATRIX == 1
%             dDOM{iprojj} = dDOM{iprojj} - BASES_RVE.DISPLACEMENTS.U*(PRB\(MtimesU'*dDOM{iprojj})) ;
%
%             %                         if ~isempty(DISPLACEMENTS_INTERFACES_BOUNDARY.MODES)
%             %                               DISPLACEMENTS_INTERFACES_BOUNDARY.MODES = DISPLACEMENTS_INTERFACES_BOUNDARY.MODES ...
%             %                                   - BASES_RVE.DISPLACEMENTS.U*PRB*(PRB\(MtimesU'*DISPLACEMENTS_INTERFACES_BOUNDARY.MODES)) ; %(MtimesU\dDOM{iprojj}) ;
%             %                         end
%             %
%
%
%         else
%             dDOM{iprojj} = dDOM{iprojj} - BASES_RVE.DISPLACEMENTS.U*(BASES_RVE.DISPLACEMENTS.U\dDOM{iprojj}) ;
%         end
%
%         %%%%%%%%%%%
%         reactDOM{iprojj} = reactDOM{iprojj} - BASES_RVE.REACTIONS.U*(BASES_RVE.REACTIONS.U\reactDOM{iprojj}) ;
%         stressDOM{iprojj} = stressDOM{iprojj} - BASES_RVE.STRESSES.U*(BASES_RVE.STRESSES.U\stressDOM{iprojj}) ;
%
%     end
%     DATAIN.TypeMode = [' localeff'] ;
%     DATAIN.NMODES_TRUNCATE = [] ;
%     [BASES_ADD.DISPLACEMENTS,BASES_ADD.REACTIONS,BASES_ADD.STRESSES] = ...
%         DeterminationModes(DATAIN,dDOM,reactDOM,...
%         stressDOM,DATA_REFMESH,ISNOTRVE) ;
%
%     nmodesR  = size(BASES_ADD.REACTIONS.U,2) ;
%     nmodesU = size(BASES_ADD.DISPLACEMENTS.U,2) ;
%     nmodesR = min(nmodesR,nmodesU) ;
%
%     BASES.DISPLACEMENTS.U = [BASES_RVE.DISPLACEMENTS.U,  BASES_ADD.DISPLACEMENTS.U] ;
%     BASES.DISPLACEMENTS.S = [BASES_RVE.DISPLACEMENTS.S;  BASES_ADD.DISPLACEMENTS.S ;] ;
%
%
%     BASES.REACTIONS.U = [BASES_RVE.REACTIONS.U,  BASES_ADD.REACTIONS.U(:,1:nmodesR)] ;
%     BASES.REACTIONS.S = [BASES_RVE.REACTIONS.S;  BASES_ADD.REACTIONS.S(1:nmodesR)] ;
%
%     %     BASES.REACTIONS.U = [BASES_RVE.REACTIONS.U,  BASES_ADD.REACTIONS.U] ;
%     %     BASES.REACTIONS.S = [BASES_RVE.REACTIONS.S;  BASES_ADD.REACTIONS.S] ;
%
%     BASES.STRESSES.U = [BASES_RVE.STRESSES.U,  BASES_ADD.STRESSES.U] ;
%     BASES.STRESSES.S = [BASES_RVE.STRESSES.S;  BASES_ADD.STRESSES.S]  ;
%
%
% end

% Disable 11-June-2019
%BASES.DISPLACEMENTS_INTERFACES_BOUNDARY = DISPLACEMENTS_INTERFACES_BOUNDARY ;

end



