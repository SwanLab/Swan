clc
clear all
% Modal analysis domain-wise. See DomainDecom_SVD.pdf
addpath('DATA_input') ;
addpath('FE_CODE') ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT_DATAFILE  ={'DATA_BEAMporousBEND'} ;  {'DATA_BEAMsolidFINE'} ; ; {'DATA_BEAM_Hexa2D_D6_bend'} ; {'DATA_BEAMporousBEND'};;{'DATA_BEAMsolidENA_otherEND'} ; {'DATA_BEAMsolidENA','DATA_BEAMsolidENA_otherEND'} ; {'DATA_BEAM_Hexa2D_D6_bend'} ;  'DATA_BEAM_Hexa2D_bend' ;  % MACRO-STRUCTURE
INPUT_DATA_RVE = [];
NROWS = [4] ;
NCOLS = [20] ;
ROWS =  { };
COLS =  {} ;
DATA.TypeUnitCell = 'HEXAG_2D_SQUARE';    % Type of unit cell
NMODES_SHOW = [] ;
SVD_FACES = 0 ;
DATA.PLOT_RECONSTRUCTED_DISPLACEMENTS = 1 ;  % Only valid for one project
DATA.NMODES_TRUNCATE = [2] ;   % Number of modes for truncation
%
% INCLUDE_MODES_RVE = 0 ;
COORrve = {} ;
CNrve = {}  ;
dRVE = {}  ;
NODESrve = {}  ;
NODESfaces = {}  ;
RIGID_BODY_MOTIONglo = {} ;
CoordinatesChangeGLO = {} ;

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
    % outputFErve = ['DATAWS/',INPUT_DATA_RVE,'_WS.mat'] ;
    
    %nameDISP =  ['DATAWS/DISP_',INPUT_DATA_RVE,'.mat'] ;
    
    %%%%%%%%%%%%%%%%%%%%%
    %     if  INCLUDE_MODES_RVE == 1
    %         load(outputFErve,'COOR','CN') ;
    %         COOR_cell  = COOR ;
    %         CN_cell = CN ;
    %     end
    load(outputFE,'MaterialType','COOR','CN','d','TypeElement','posgp') ;
    % Number of materials/bodies
    IND_MAT =unique(MaterialType) ;
    nRVE = length(IND_MAT) ;
    
    
    % Store the coordinates, connectivitites and displacememts of each cell in
    % a cell array
    ndim = size(COOR,2) ;
    d = reshape(d,ndim,[])' ; % Displacement vector
    
    % Tolerance for matching nodes
    TOL = ChooseTolerance(CN,COOR) ;
    
    
    
    for i = 1:nRVE
        IndElements =  find(IND_MAT(i) == MaterialType) ;
        CNrve{end+1} =  CN(IndElements,:) ; % Connectivities
        NODESrve{end+1} =unique( CNrve{end}(:)) ;  % Nodes
        % Coordinates
        COORabs = COOR(NODESrve{end},:) ;
        
        xmin = min(COORabs(:,1)) ; xmin = xmin(1) ;
        xmax = max(COORabs(:,1)) ; xmax = xmax(1) ;
        ymin = min(COORabs(:,2)) ; ymin = ymin(1) ;
        ymax = max(COORabs(:,2)) ; ymax = ymax(1) ;
        zmin = [] ;   zmax = [] ;
        dRVE{end+1} = d( NODESrve{end},:) ;
        
        
        % Boundaries of each cell
        switch  DATA.TypeUnitCell
            case 'HEXAG_2D_SQUARE'
                [NODESfaces{end+1} NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORabs,TOL,xmax,xmin,ymax,ymin,zmax,zmin,1:size(COORabs,1)) ;
                % REFERENCE POINT  (TO MEASURE COORDINATES, AS WELL AS DISPLACEMENTS)
                % BOTTOM NODES, leftmost
                lbott = NODESfaces{end}{2} ;
                [xminREF NODEREF] = min(COORabs(lbott,1)) ;
                NODEREF = lbott(NODEREF) ;
                % Reference point 2  (top)
                lbott = NODESfaces{end}{4 } ;
                [xminREF NODEREF2] = min(COORabs(lbott,1)) ;
                NODEREF2 = lbott(NODEREF2) ;
                
            otherwise
                error('Option not implemented')
        end
        % Displacements
        
        % Now coordinates and displacements are referred to the reference point
        nnode = size(COORabs,1) ;
        RIGID_BODY_MOTION = repmat(dRVE{end}(NODEREF,:),nnode,1) ;
        CoordinatesChange = repmat(COORabs(NODEREF,:),nnode,1);  % Change of coordinates
        COORrve{end+1} = COORabs - CoordinatesChange ;
        %  dRVE{end} =  dRVE{end} -RIGID_BODY_MOTION  ;  % Avoid translations
        % Purging rotations
        %  dispHORIZONTAL2 = dRVE{end}(NODEREF2,1) ;  % Horizontal displacement 2nd reference point
        %   distanceBETWEEN_refP = (COORrve{end}(NODEREF2,2) - COORrve{end}(NODEREF,2))  ;
        %   ANGLE_ROTATION = dispHORIZONTAL2/distanceBETWEEN_refP;
        
        %    disp_HORIZONTAL = ANGLE_ROTATION*COORrve{end}(:,2) ;
        
        %   dRVE{end}(:,1) =  dRVE{end}(:,1) - disp_HORIZONTAL ;  % Avoid rotations
        
        %  RIGID_BODY_MOTION(:,1) =    RIGID_BODY_MOTION(:,1) + disp_HORIZONTAL ;
        
        %   RIGID_BODY_MOTION = reshape(RIGID_BODY_MOTION',[],1) ;
        
        %  RIGID_BODY_MOTIONglo{end+1} = RIGID_BODY_MOTION ;
        CoordinatesChangeGLO{end+1} = CoordinatesChange ;
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%   unit cell modes
    %     if  INCLUDE_MODES_RVE == 1
    %         load(nameDISP,'DISPLACEMENTS') ;
    %         d = DISPLACEMENTS ;
    %         for i = 1:size(d,2)
    %             CNrve{end+1} =  CN_cell ; % Connectivities
    %             NODESrve{end+1} = 1:size(COOR_cell,1)' ;
    %             % Coordinates
    %             COORabs = COOR_cell ;
    %
    %             xmin = min(COORabs(:,1)) ; xmin = xmin(1) ;
    %             xmax = max(COORabs(:,1)) ; xmax = xmax(1) ;
    %             ymin = min(COORabs(:,2)) ; ymin = ymin(1) ;
    %             ymax = max(COORabs(:,2)) ; ymax = ymax(1) ;
    %             zmin = [] ;   zmax = [] ;
    %             d = reshape(DISPLACEMENTS(:,i),ndim,[])' ;
    %
    %             dRVE{end+1} = d;
    %
    %
    %             % Boundaries of each cell
    %             switch  DATA.TypeUnitCell
    %                 case 'HEXAG_2D_SQUARE'
    %                     [NODESfaces{end+1} NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORabs,TOL,xmax,xmin,ymax,ymin,zmax,zmin,1:size(COORabs,1)) ;
    %                     % REFERENCE POINT  (TO MEASURE COORDINATES, AS WELL AS DISPLACEMENTS)
    %                     % BOTTOM NODES, leftmost
    %                     lbott = NODESfaces{end}{2} ;
    %                     [xminREF NODEREF] = min(COORabs(lbott,1)) ;
    %                     NODEREF = lbott(NODEREF) ;
    %                     % Reference point 2  (top)
    %                     lbott = NODESfaces{end}{4 } ;
    %                     [xminREF NODEREF2] = min(COORabs(lbott,1)) ;
    %                     NODEREF2 = lbott(NODEREF2) ;
    %
    %                 otherwise
    %                     error('Option not implemented')
    %             end
    %             % Displacements
    %
    %             % Now coordinates and displacements are referred to the reference point
    %             nnode = size(COORabs,1) ;
    %             COORrve{end+1} = COORabs - repmat(COORabs(NODEREF,:),nnode,1) ;
    %             dRVE{end} =  dRVE{end} - repmat(dRVE{end}(NODEREF,:),nnode,1) ;  % Avoid translations
    %             % Purging rotations
    %             dispHORIZONTAL2 = dRVE{end}(NODEREF2,1) ;  % Horizontal displacement 2nd reference point
    %             distanceBETWEEN_refP = (COORrve{end}(NODEREF2,2) - COORrve{end}(NODEREF,2))  ;
    %             ANGLE_ROTATION = dispHORIZONTAL2/distanceBETWEEN_refP;
    %
    %             disp_HORIZONTAL = ANGLE_ROTATION*COORrve{end}(:,2) ;
    %
    %             dRVE{end}(:,1) =  dRVE{end}(:,1) - disp_HORIZONTAL ;  % Avoid rotations
    %
    %
    %         end
    %     end
    
end

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
for i = 1:length(COORrve)
    if i ~= refMESH
        IDX = knnsearch(COORrve{i},COORref) ;
        COORrve{i} =  COORrve{i}(IDX,:) ;
        CoordinatesChange=  CoordinatesChangeGLO{i}(IDX,:) ;
        CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChange)),1 );
        COORnew{i} =COORrefORD + CoordinatesChange ;
        
        %     RIGID_BODY_MOTIONglo{i} =  RIGID_BODY_MOTIONglo{i}(IDX,:) ;
        %     RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
        dRVE{i} =  dRVE{i}(IDX,:)  ;
        dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
        for j = 1:length(NODESfaces{i})
            NODESfacesALL = zeros(size(NODESrve{i})) ;
            NODESfacesALL(NODESfaces{i}{j}) = 1 ;
            
            NODESfacesALL = NODESfacesALL(IDX) ;
            NODESfaces{i}{j}  =find(NODESfacesALL == 1) ;
            
        end
        
        NODESrve{i} =  NODESrve{i}(IDX) ;
        
    else
        dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
        %      RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
        CoordinatesChange=  CoordinatesChangeGLO{i}  ;
        CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChangeGLO{i})),1 );
        COORnew{i} =COORrefORD + CoordinatesChange ;
    end
end








%% SVD of displacement field
addpath('SVDlibrary')

%%%% MODES UNIT CELL
%if ~isempty(INDICES_CELL)
%    SNAP = cell2mat(dRVE(INDICES_CELL)) ;
%else
SNAP = cell2mat(dRVE) ;


%%%%%%%%%%%%%%%%%%%%%%%
if ndim ==2 ; nRB = 3;  else ; error('Not implemented') ; end ;
NODESref = NODESrve{refMESH}  ;
nnodeRVE = length(NODESref) ;
MODErb = zeros(size(SNAP,1),nRB) ;
MODErb(1:ndim:end,1) = 1/sqrt(nnodeRVE) ;
MODErb(2:ndim:end,2) = 1/sqrt(nnodeRVE) ; ;
COORref = COOR(NODESref,:) ;
xGEO = sum(COORref(:,1))/(nnodeRVE) ;
yGEO = sum(COORref(:,2))/(nnodeRVE) ;
COORgeo = repmat([xGEO,yGEO],nnodeRVE,1) ;
xxx = COORref(:,1)-COORgeo(:,1) ;
yyy = COORref(:,2)-COORgeo(:,2) ;
MODErb(1:ndim:end,3) = -yyy;
MODErb(2:ndim:end,3) = xxx;
MODErb(:,3) = MODErb(:,3)/norm( MODErb(:,3)) ;


%--------------------------------------
%

SNAPrg =  MODErb*(MODErb'*SNAP) ; 
SNAP = SNAP - SNAPrg;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%end

%if  INCLUDE_MODES_RVE == 1
%  %   Uhom = cell2mat(dRVE(NCELLS+1:end)) ;
%
%     %%%%%
%
%
%     %
%     errorHOMOG = 0 ;
% end
% %%%%5
%

%%%%
%
% SVD_MATLAB = 0;
%
% if SVD_MATLAB ==1
%     [U,S,V] = svd(SNAP,0) ;
%     nmodes= DATA.NMODES_TRUNCATE ;
%       if isempty(DATA.NMODES_TRUNCATE)
%           nmodes= rank(SNAP) ;
%       end
%     SNAP_recontru = U(:,1:nmodes)*(S(1:nmodes,1:nmodes)*V(:,1:nmodes)') ;
%
% else

SNAP_recontru = [] ;
[U,S,V] = SVDT(SNAP) ;
% if  INCLUDE_MODES_RVE == 1
%     U = [Uhom U ] ; % Total basis matrix
% end

SingVsq =  (S.*S) ;
nTOTAL = sqrt(sum(SingVsq) ); % + errorHOMOG^2)  ;
SingVsq = sort(SingVsq);
normEf2 = sqrt(cumsum(SingVsq)) ;
normEf2 = sort(normEf2,'descend');
hold on
lS = log10(normEf2);
% figure(1)
% hold on
% h = plot([1:length(S)-1],lS(2:end),'r');
% xlabel('MODES')
% ylabel('LOG10 SVD ERROR')
% AAA =axis ;
% if ~isempty(NMODES_SHOW)
% AAA(2)  =NMODES_SHOW ;
% axis(AAA) ;
% end

figure(2)
hold on
% firstpoint = sqrt(errorHOMOG) ;
% if  INCLUDE_MODES_RVE == 1
%     h = plot([size(Uhom,2)   (size(Uhom,2)+1):(length(S)+size(Uhom,2))],[firstpoint; normEf2]/nTOTAL*100,'r*');
%     xlabel('MODES')
%     ylabel('SVD ERROR  (%) (additional homog.)')
%
% else
h = plot([ 1:length(S)-1],[ normEf2(2:end)]/nTOTAL*100,'r-*');
xlabel('MODES')
ylabel('SVD ERROR  (%) ')
AAA =axis ;
if ~isempty(NMODES_SHOW)
    AAA(2)  =NMODES_SHOW ;
    axis(AAA) ;
end

%end
% end

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
        METHOD  = 0;
        if METHOD == 1
            
            Vtrun = bsxfun(@times,Vtrun',Strun)' ;
            SNAP_recontru = Utrun*Vtrun' ;   % Reconstructed
            
        else
            SNAP_recontru = Utrun*(diag(Strun)*Vtrun') ;
        end
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
    %    RIGID_BODY_MOTIONglo =cell2mat(RIGID_BODY_MOTIONglo) ;
    
    % Just out of curiosity: what is the rank of RIGID_BODY_MOTIONglo ?
    %    [MODESrb,Srb,Vrb] = SVDT(RIGID_BODY_MOTIONglo,0) ;
    %  figure(85)
    %  hold on
    %  xlabel('Modes RIGID BODY')
    
    %     SingVsq =  (Srb.*Srb) ;
    %     nTOTAL = sqrt(sum(SingVsq) ); % + errorHOMOG^2)  ;
    %     SingVsq = sort(SingVsq);
    %     normEf2 = sqrt(cumsum(SingVsq)) ;
    %     normEf2 = sort(normEf2,'descend');
    %     h = plot([ 1:length(Srb)-1],[ normEf2(2:end)]/nTOTAL*100,'r-*');
    %     xlabel('MODES')
    %     ylabel('SVD ERROR RIGID BODY MODES  (%) ')
    %     AAA =axis ;
    %     AAA(2)  =NMODES_SHOW ;
    %     axis(AAA) ;
    %
    %     AAA(2)  =NMODES_SHOW ;
    %     axis(AAA) ;
    %
    %
    
    dRECONS = SNAP_recontru + SNAPrg ;
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
end


%%%%%%%%%%%%%
% SVD displacement faces parallel to x = 0
% -----------------------------------------
%refMESH
%NODESfaces
if SVD_FACES == 1
    switch DATA.TypeUnitCell
        case 'HEXAG_2D_SQUARE'
            NODES= NODESfaces{refMESH} ;
            NODESref = NODESrve{refMESH} ;
            % Dispalcements Faces parallel to x-axix
            SNAPx = [U(NODES{1},1:NMODES_SHOW) U(NODES{3},1:NMODES_SHOW)] ;
            % Dispalcements Faces parallel to x-axix
            
            [Ux,S,V] = SVDT(SNAPx) ;
            % if  INCLUDE_MODES_RVE == 1
            %     U = [Uhom U ] ; % Total basis matrix
            % end
            
            SingVsq =  (S.*S) ;
            nTOTAL = sqrt(sum(SingVsq) ); % + errorHOMOG^2)  ;
            SingVsq = sort(SingVsq);
            normEf2 = sqrt(cumsum(SingVsq)) ;
            normEf2 = sort(normEf2,'descend');
            figure(3)
            hold on
            
            h0 = plot([ 1:length(S)-1],[ normEf2(2:end)]/nTOTAL*100,'g-*');
            xlabel('MODES')
            ylabel('SVD ERROR -xFACE (%) ')
            AAA =axis ;
            AAA(2)  =NMODES_SHOW ;
            axis(AAA) ;
            % end
            
            SNAPy = [U(NODES{2},1:NMODES_SHOW) U(NODES{4},1:NMODES_SHOW)] ;
            [Ux,S,V] = SVDT(SNAPy) ;
            % if  INCLUDE_MODES_RVE == 1
            %     U = [Uhom U ] ; % Total basis matrix
            % end
            
            SingVsq =  (S.*S) ;
            nTOTAL = sqrt(sum(SingVsq) ); % + errorHOMOG^2)  ;
            SingVsq = sort(SingVsq);
            normEf2 = sqrt(cumsum(SingVsq)) ;
            normEf2 = sort(normEf2,'descend');
            figure(3)
            hold on
            
            h1 = plot([ 1:length(S)-1],[ normEf2(2:end)]/nTOTAL*100,'r-*');
            xlabel('MODES')
            ylabel('SVD ERROR -yFACE (%) ')
            AAA =axis ;
            AAA(2)  =NMODES_SHOW ;
            axis(AAA) ;
            % end
            
            
            legend([h0,h1],{'xFACE','yFACE'})
    end
    
end

%%%%%%%




%%%% PLOTTING MODES
CNref = CNrve{refMESH} ;
NODESref = NODESrve{refMESH}  ;
DOFl = [] ;
DATA.NODES = NODESref' ;
MODES= U ;
GidPostProcessModes(COOR,CNref,TypeElement,MODES,posgp,NAME_MODES,DATA,DOFl);

%%%% PLOTTING MODES
%  
% MODES= MODErb ;
% NAME_MODES = ['RIGID_',NAME_MODES];
% GidPostProcessModes(COOR,CNref,TypeElement,MODES,posgp,NAME_MODES,DATA,DOFl);
