clc
clear all
% Modal analysis domain-wise,    See DomainDecom_SVD.pdf
addpath('DATA_input') ;
addpath('FE_CODE') ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT_DATAFILE  ='DATA_BEAM3Drep'; 'DATAStent1' ;  'DATA_BEAMhexaIND'; 'DATA_BEAMhexa'; 'DATA_BEAMthinw';  'DATA_BEAMsolidEMP30'  ;'DATA_BEAMslice'   ;; 'DATA_BEAMsolid' ;{'DATA_BEAMporousBEND'} ; {'DATA_BEAMporousBEND'}; {'DATA_BEAM_Hexa2D_D6_bend'} ; {'DATA_BEAMporousBEND'};;{'DATA_BEAMsolidENA_otherEND'} ; {'DATA_BEAMsolidENA','DATA_BEAMsolidENA_otherEND'} ; {'DATA_BEAM_Hexa2D_D6_bend'} ;  'DATA_BEAM_Hexa2D_bend' ;  % MACRO-STRUCTURE
DATA.NMODES_TRUNCATE = [10] ;   % Number of modes for truncation
DATA.NMODES_TRUNCATE_RIGID_BODY = [] ;   % Number of modes for truncation of rigid body
DATA.METHOD_PURGE_ROTATIONS = 'SVD_RIGID_DEF' ;
DATA.NAMEWS = ['DATAWS/',INPUT_DATAFILE,'_OFFLINE','.mat'] ; % Binary file where information is saved
DATA.METHOD_COMPUTE_REACTION_BASIS = 'SVD' ;'RIGHT_SPACE' ;    ;'SVD_PROJECT' ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INPUT_DATA_RVE = [];
NROWS = [4] ;
NCOLS = [20] ;
ROWS =  { };
COLS =  {} ;
DATA.TypeUnitCell = 'HEXAG_2D_SQUARE';    % Type of unit cell
NMODES_SHOW = 20 ;
SVD_FACES = 0 ;
DATA.PLOT_RECONSTRUCTED_DISPLACEMENTS = 1 ;  % Only valid for one project

%
% INCLUDE_MODES_RVE = 0 ;
COORrve = {} ;
CNrve = {}  ;
dRVE = {}  ;
NODESrve = {}  ;
NODESfaces = {}  ;
RIGID_BODY_MOTIONglo = {} ;
CoordinatesChangeGLO = {} ;
reactDOM = {} ; % REactions of each domain


if ~iscell(INPUT_DATAFILE)
    INPUT_DATAFILE = {INPUT_DATAFILE}  ;
end
IndElements = {} ;

for iproject = 1:length(INPUT_DATAFILE)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    NAME_MODES =[INPUT_DATAFILE{iproject},'.msh'];
    %%%%
    if ~isempty(ROWS)
        INDICES_CELL  = [];
        for irows = 1:length(ROWS{iproject})
            INDICES_CELLloc = (ROWS{iproject}(irows)-1)*NCOLS(iproject) +COLS{iproject};
            INDICES_CELL = [INDICES_CELL INDICES_CELLloc] ;
        end
    else
        INDICES_CELL =[] ;
    end
    
    
    outputFE = ['DATAWS/',INPUT_DATAFILE{iproject},'_WS.mat'] ;
    
    load(outputFE,'MaterialType','COOR','CN','d','TypeElement','posgp','Nst','wSTs_RHS','fNOD','posgp_RHS',...
        'CNb','TypeElementB','Fpnt','Tnod','CONNECTb','stressGLO','wSTs','Bst') ;
    % Number of materials/bodies
    IND_MAT =unique(MaterialType) ;
    nRVE = length(IND_MAT) ;
    
    
    % Store the coordinates, connectivitites and displacememts of each cell in
    % a cell array
    ndim = size(COOR,2) ;
    d = reshape(d,ndim,[])' ; % Displacement vector
    
    % Tolerance for matching nodes
    TOL = ChooseTolerance(CN,COOR) ;
    
    
    % Loop over number of DOMAINS
    for i = 1:nRVE
        IndElements{i} =  find(IND_MAT(i) == MaterialType) ;   % Elements of domain i
        CNrve{i} =  CN(IndElements{i},:) ; % Connectivities of Domain i
        NODESrve{i} =unique( CNrve{i}(:)) ;  % Nodes of Domain i
        COORabs = COOR(NODESrve{i},:) ;  % Coordinates of DOmain i
        
        dRVE{i} = d( NODESrve{i},:) ;  % Displacements of domain i
        
        % Purging rotations
        switch DATA.METHOD_PURGE_ROTATIONS
            case 'SVD_RIGID_DEF'
                [NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,COORrve] ...
                    = PurgingRotations(i,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
                    COORrve);
            case 'ANALYTICAL'
                % No longer used
                error('This option is no longer available')
                [NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,COORrve] ...
                    = PurgingRotationsANALY(i,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
                    COORrve,NODESrve,CNrve);
        end
        
        %%% Computing self-equilibrated reactions for each domain
        % -------------------------------------------------------
        [reactDOMloc] = ReactionsRVE(i,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
            COORrve,Nst,wSTs_RHS,fNOD,NODESrve,CNrve,posgp,IndElements{i},posgp_RHS,...
            CNb,TypeElementB,Fpnt,Tnod,CONNECTb,COOR,stressGLO,wSTs,Bst,TypeElement);
        reactDOM{i} =  reactDOMloc  ;
        
        
    end
    
    
    
    
    
end

dbstop('116')

NCELLS = length(dRVE) ;

%%% REFERENCE MESH. We choose as reference mesh that with less number of
%%% elements
[nelems,ncols] = cellfun(@size,CNrve) ;
[nelemMAX, refMESH] = min(nelems) ;
%% COORDINATES OF REFERENCE MESH
COORref = COORrve{refMESH} ;
COORrefORD = reshape(COORref',prod(size(COORref)),1) ;
% Now we re-order
COORnew = {} ;
ElementsREF = IndElements{refMESH} ;
for i = 1:length(COORrve)
    if i ~= refMESH
        IDX = knnsearch(COORrve{i},COORref) ;
        COORrve{i} =  COORrve{i}(IDX,:) ;
        CoordinatesChange=  CoordinatesChangeGLO{i}(IDX,:) ;
        CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChange)),1 );
        COORnew{i} =COORrefORD + CoordinatesChange ;
        
        RIGID_BODY_MOTIONglo{i} =  RIGID_BODY_MOTIONglo{i}(IDX,:) ;
        RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
        dRVE{i} =  dRVE{i}(IDX,:)  ;
        dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
        IDXdof = small2large(IDX,ndim) ;
        reactDOM{i} = reactDOM{i}(IDXdof) ;
        for j = 1:length(NODESfaces{i})
            NODESfacesALL = zeros(size(NODESrve{i})) ;
            NODESfacesALL(NODESfaces{i}{j}) = 1 ;
            
            NODESfacesALL = NODESfacesALL(IDX) ;
            NODESfaces{i}{j}  =find(NODESfacesALL == 1) ;
            
        end
        
        NODESrve{i} =  NODESrve{i}(IDX) ;
        
    else
        dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
        RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
        CoordinatesChange=  CoordinatesChangeGLO{i}  ;
        CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChangeGLO{i})),1 );
        COORnew{i} =COORrefORD + CoordinatesChange ;
    end
end

%% SVD of displacement field
addpath('SVDlibrary')
nfigure = 1;
LEGENDG = 'SVD  ERROR (%)' ;
COLOR = 'r' ;
[U,S,V,h1] = SVD_and_error(dRVE,nfigure,LEGENDG,NMODES_SHOW,COLOR ) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% SVD of reaction forces
nfigure = 1;
LEGENDG = 'SVD  ERROR (%)' ;
COLOR = 'b' ;

switch DATA.METHOD_COMPUTE_REACTION_BASIS
    case  'SVD_PROJECT'
        hold on
        legend([h1 ],{'DISPLACEMENTS'})
        BasisUdef = U(:,1:DATA.NMODES_TRUNCATE) ; % = [2] ;
        
        NODES_FACES = NODESfaces{refMESH} ;
        f1 = NODES_FACES{1} ;
        f2 = NODES_FACES{3} ;
        f1 = small2large(f1,ndim) ;
        f2 = small2large(f2,ndim) ;
        f = [f1; f2] ;
        reactDOM  =  cell2mat(reactDOM) ;
        reactDOM = reactDOM(f,:) ;
        reactDOM = BasisUdef(f,:)*(BasisUdef(f,:)\reactDOM) ;
        [UR,SR,VR] = SVDT(reactDOM,0) ;
        
        BasisRdef = zeros(size(BasisUdef)) ;
        BasisRdef(f,:) = UR ;
        figure(2)
        hold on
        xlabel('Reaction Modes')
        ylabel('Singular Values')
        plot(SR)
        
    case 'SVD'
        [BasisRdef,Sr,Vr,h2] = SVD_and_error(reactDOM,nfigure,LEGENDG,NMODES_SHOW,COLOR ) ;
        % BasisRdef = BasisRdef(:,1:DATA.NMODES_TRUNCATE);
        %   BasisUdef = U(:,1:DATA.NMODES_TRUNCATE) ; % = [2] ;
        hold on
        try
            legend([h1 h2],{'DISPLACEMENTS','REACTIONS'})
        end
        
    case 'RIGHT_SPACE'
        BasisRdef = zeros(size(U)) ;
        reactDOM = cell2mat(reactDOM) ;
        for imode =1:size(U,2)
            BasisRdef(:,imode) = reactDOM*V(:,imode)/(norm(reactDOM*V(:,imode))) ;
        end
        
end


save(DATA.NAMEWS,'BasisRdef','ElementsREF')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOTING RECONSTRUCTED DISPLACEMENTS
% ------------------------------------
if  DATA.PLOT_RECONSTRUCTED_DISPLACEMENTS == 1 ;  % Only valid for one project
    if length(INPUT_DATAFILE) >1
        error('Option not implemented yet')
    end
    %NMODES_TRUNCATE
    % Truncated displacement field
    %if isempty(SNAP_recontru)
    if isempty(DATA.NMODES_TRUNCATE)
        SNAP_recontru = SNAP ;
    else
        
        
        Utrun = U(:,1:DATA.NMODES_TRUNCATE) ;
        Strun = S(1:DATA.NMODES_TRUNCATE) ;
        Vtrun = V(:,1:DATA.NMODES_TRUNCATE) ;
        
        
        save(DATA.NAMEWS,'U','S','V','-append');
        
        
        %         METHOD  = 0;
        %         if METHOD == 1
        
        Vtrun = bsxfun(@times,Vtrun',Strun)' ;
        SNAP_recontru = Utrun*Vtrun' ;   % Reconstructed
        
        %         else
        %             SNAP_recontru = Utrun*(diag(Strun)*Vtrun') ;
        %         end
        %else
        
    end
    % end
    % Next we sum the rigid body motions
    %
    %  Each column of SNAP_reconstru represents the displacement field of
    %  the corresponding unit cell. Accordingly, all we have to do is to
    %  place all displacements in a single vector
    %  dRECONS =SNAP_recontru(:) ; % This is our new "displacement vector" associated to deformation.
    % We have to sum up the rigid body motions.
    RIGID_BODY_MOTIONglo =cell2mat(RIGID_BODY_MOTIONglo) ;
    
    % Just out of curiosity: what is the rank of RIGID_BODY_MOTIONglo ?
    [MODESrb,Srb,Vrb] = SVDT(RIGID_BODY_MOTIONglo,0) ;
    figure(85)
    hold on
    xlabel('Modes RIGID BODY')
    
    SingVsq =  (Srb.*Srb) ;
    nTOTAL = sqrt(sum(SingVsq) ); % + errorHOMOG^2)  ;
    SingVsq = sort(SingVsq);
    normEf2 = sqrt(cumsum(SingVsq)) ;
    normEf2 = sort(normEf2,'descend');
    h = plot([ 1:length(Srb)-1],[ normEf2(2:end)]/nTOTAL*100,'r-*');
    xlabel('MODES')
    ylabel('SVD ERROR RIGID BODY MODES  (%) ')
    AAA =axis ;
    AAA(2)  =NMODES_SHOW ;
    axis(AAA) ;
    
    AAA(2)  =NMODES_SHOW ;
    axis(AAA) ;
    
    
    if isempty(DATA.NMODES_TRUNCATE_RIGID_BODY )
        SNAPRB_RECONS =    RIGID_BODY_MOTIONglo  ;
    else
        Utrun = MODESrb(:,1:DATA.NMODES_TRUNCATE_RIGID_BODY) ;
        Strun = Srb(1:DATA.NMODES_TRUNCATE_RIGID_BODY) ;
        Vtrun = Vrb(:,1:DATA.NMODES_TRUNCATE_RIGID_BODY) ;
        %         METHOD  = 0;
        %         if METHOD == 1
        
        Vtrun = bsxfun(@times,Vtrun',Strun)' ;
        SNAPRB_RECONS = Utrun*Vtrun' ;   % Reconstructed
    end
    
    
    
    dRECONS = SNAP_recontru + SNAPRB_RECONS ;
    % And what about our coordinate matrix  ? Information is contained in
    % COORnew.
    COORnew = cell2mat(COORnew) ;
    %The difficulty lies in the CONNECTIVITY MATRIX. How to do it
    % ? The required information is in
    % CNref = CNrve{refMESH} ;
    % NODESref = NODESrve{refMESH}  ;
    % How to handle it ? Intuitively, we have to make several "translated"
    % copies of  CNref matrix.
    % We begin by renumbering CNref
    CNref = CNrve{refMESH} ;
    CNrefNEW = zeros(size(CNref)) ;
    NODESref = NODESrve{refMESH}  ;
    for inode = 1:length(NODESref)
        nodeLOC = NODESref(inode) ;
        INDnodes = find(CNref==nodeLOC) ;
        CNrefNEW(INDnodes) = inode ;
    end
    % Next we  reorder all matrices by placing as first mesh the one
    % corresponding to the reference element.
    REMAINING_meshes = setdiff(1:size(COORnew,2),refMESH) ;
    COORnew = COORnew(:,[refMESH REMAINING_meshes]) ;
    dRECONS =  dRECONS(:,[refMESH REMAINING_meshes]) ;
    %-------------------------------------------------------------------
    % CONSTRUCTING THE MATRIX OF COORDINATES, CONNECTIVITIES AND MATERIAL
    % LIBRARY
    % -------------------------------------------------------------------
    COORrecons = [] ;
    CNrecons = []  ;
    Materials = [] ;
    
    % Connectivities
    for e = 1:size(COORnew,2)
        COORloc = reshape(COORnew(:,e),ndim,[])'  ;
        COORrecons = [COORrecons; COORloc] ;
        CNloc = CNrefNEW +(e-1)*size(COORloc,1) ;
        CNrecons = [CNrecons ; CNloc] ;
        Materials =[Materials; e*ones(size(CNloc,1),1)]  ;
        
    end
    % Material library
    
    NAME_BASE_GIDfiles = [INPUT_DATAFILE{1},'_RECONS_nmodes','_',num2str(DATA.NMODES_TRUNCATE)] ;
    NameFileMeshss = [] ;
    dRECONS = dRECONS(:) ;
    strainGLOgid = [] ;
    stressGLOgid = [] ;
    React = [] ;
    GidPostProcess(COORrecons,CNrecons,TypeElement,dRECONS,strainGLOgid, ...
        stressGLOgid,  React,NAME_BASE_GIDfiles,posgp,NameFileMeshss,Materials,DATA);
    MaterialsRECONS = Materials ;
    %  save(DATA.NAMEWS,'COORrecons','CNrecons','MaterialsRECONS','-append')
    
end





%%%% PLOTTING MODES
CNref = CNrve{refMESH} ;
NODESref = NODESrve{refMESH}  ;
DOFl = [] ;
DATA.NODES = NODESref' ;

MODES= U ;
disp('Displacements')

GidPostProcessModes(COOR,CNref,TypeElement,MODES,posgp,NAME_MODES,DATA,DOFl);

save(DATA.NAMEWS,'-append','CNref','NODESref','COOR','TypeElement','posgp');

%%%% PLOTTING MODES
disp('Reactions')
MODES= BasisRdef;
NAME_MODES = ['REACTIONS',NAME_MODES] ;
GidPostProcessModes(COOR,CNref,TypeElement,MODES,posgp,NAME_MODES,DATA,DOFl);
