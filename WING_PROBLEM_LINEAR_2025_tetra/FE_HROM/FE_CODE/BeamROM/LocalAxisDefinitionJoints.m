function [ROTATIONS,INFOLINES,SUCCESSFUL_MATCH]  = LocalAxisDefinitionJoints(MESH1D,MESH3D,DATAIN)

% ROTATION MATRICES for joints
% -----------------
% ------------------------
% Loop over entities loop
% --------------------------
if nargin == 0
    load('tmp1.mat')
end
ndim = size(MESH1D.COOR,2);
INFOLINES = MESH1D.INFOLINES ; %
ROTATIONS = MESH1D.ROTATIONS ;
SUCCESSFUL_MATCH = cell(length(MESH1D.PROP),1) ;
for imat1D  = 1: length(MESH1D.PROP)  % Loop over types of entities (both beams and joints )
    switch MESH1D.PROP(imat1D).TYPE ;  % Check whether it is a joint or a beam
        case 'JOINT'
            % It is a group of  joints
            itypeJOINT = MESH1D.PROP(imat1D).INDEX_ENTITY ;  %  index that identifies the 3D mesh
            
            if ~isempty(MESH1D.INFOLINES.NODES{imat1D})
                % Joints with just 2 ends
                END_NODES_GLO = {} ; END_ELEMENTS_GLO = {} ;
                % Loop over joints of this group
                for ijoints =1:length(MESH1D.INFOLINES.NODES{imat1D})
                    END_NODES_GLO{ijoints} = [MESH1D.INFOLINES.NODES{imat1D}{ijoints}([1,end])   ] ;
                    END_ELEMENTS_GLO{ijoints} = [MESH1D.INFOLINES.ELEMENTS{imat1D}{ijoints}([1,end]) ] ;
                    
                    %    END_ELEMENTS_GLO{ijoints} = MESH1D.INFOLINES.ELEMENTS{imat1D}([1,end]) ;
                end
            else
                % Joints with more than 2 ends
                END_NODES_GLO = MESH1D.INFOLINES.END_NODES{imat1D} ;
                END_ELEMENTS_GLO= MESH1D.INFOLINES.END_ELEMENTS{imat1D} ;
            end
            END_NODES_GLO_NEW = cell(size(END_NODES_GLO)) ;   % Goal: Recalculate the order of these nodes
            END_ELEMENTS_GLO_NEW = cell(size(END_ELEMENTS_GLO)) ;
            % ----------------------------------
            % Reference configuration JOINT
            % -----------------------------------
            CENTROIDS = MESH3D.JOINTS(itypeJOINT).DATA3D.CENTROIDS ;  % Centrois of all faces 
            NODES_FACES_REFERENCE = MESH3D.JOINTS(itypeJOINT).DATA3D.NODES_FACES ; % Nodes of all faces
            NCONNECTING_FACES = length(size(END_NODES_GLO{1})) ; 
            COOR_FACES_REFERENCE = cell(NCONNECTING_FACES,1) ; % Coordinates of connecting faces 
            for iface = 1:length(COOR_FACES_REFERENCE)
                COOR_FACES_REFERENCE{iface} = MESH3D.JOINTS(itypeJOINT).DATA3D.COOR(NODES_FACES_REFERENCE{iface},:)' ;
            end
            CENTROIDS = CENTROIDS(:,1:NCONNECTING_FACES) ; 
            % Loop over "joints" of the  group under consideration
            % -----------------------
            SUCCESSFUL_MATCH{imat1D}  = zeros(1,length(END_NODES_GLO)) ;
            for ijoint = 1:length(END_NODES_GLO)
                % 1D mesh
                % -------
                END_NODES = END_NODES_GLO{ijoint} ;    % End nodes (1D mesh)
                END_ELEMENTS = END_ELEMENTS_GLO{ijoint} ;    % End elements (1D mesh)
                %                  END_NODES = MESH1D.INFOLINES.END_NODES{ientity}
                
                
                %                   if isempty(END_NODES)
                %                     % Two-faces joint
                %                     % ---------------
                %                     END_NODES = MESH1D.INFOLINES.NODES{ientity}{SELECTED_JOINT}([1,end]) ;
                %                     END_ELEMENTS = MESH1D.INFOLINES.ELEMENTS{ientity}{SELECTED_JOINT}([1,end]) ;
                %                 end
                
                
                
                COOR_glo  = MESH1D.COOR(END_NODES,:)' ; % Coordinates end nodes
                % -------------------------------------------
                % Coordinate nodes faces of the connecting slices  (neighboring slices)
                % ----------------------------------------
                CONNECTslices = ConnectingSlices_JOINT(END_NODES,END_ELEMENTS,MESH1D,MESH3D) ;
%  For instance                CONNECTslices = 
% 
%               COOR: {[3x256 double]  [3x256 double]}  % Coordinate face
%               nodes
%              NODES: {[256x1 double]  [256x1 double]}
%        TRANSLATION: {[3x1 double]  [3x1 double]}  .... Translation that
%        has to undergo the slice to be placed in the position dicated by
%        the 1D MESH
%           ROTATION: {[3x3 double]  [3x3 double]}
%     indexSLICE_GLO: [1 1]
                
                
                % All we have to do is to find the rotation and translation that makes that coordinates
                %COOR_FACES_REFERENCE match CONNECTslices.COOR.
                TOL = MESH3D.JOINTS(itypeJOINT).DATA3D.TOLERANCE_DISTANCE_NODES;
                % Rotation is defined with respect to the centroid of face
                % F1 (which face is F1 is determined by the user via GID )
                [FACES_ORDER, ROTATION_FACE1,SUCCESSFUL_MATCH_LOC ] =  ...
                    MatchFacesJointSlice(COOR_glo,CONNECTslices.COOR,CENTROIDS,COOR_FACES_REFERENCE,TOL) ;
                SUCCESSFUL_MATCH{imat1D}(ijoint)  = SUCCESSFUL_MATCH_LOC ;
                % New order of list "END_NODES" and "END_ELEMENTS"
                END_NODES_GLO_NEW{ijoint} = END_NODES(FACES_ORDER) ;
                END_ELEMENTS_GLO_NEW{ijoint} = END_ELEMENTS(FACES_ORDER) ;
                %% STORE ROTATION MATRICES IN THEIR CORRESPONDING entries
                % Recall that
                %   these are the rotations between the configuration of the
                % joint in the structure, and the configuration in which
                % the joint was defined (meshed).
                % For joints, this is defined with respect to the first
                % local node. The remaining 1D elements have no rotation
                % matrix assigned
                ielem = END_ELEMENTS_GLO_NEW{ijoint}(1) ;
                ifin = ielem*ndim ; iini = ifin-ndim+1;
                ROTATIONS(:,iini:ifin) = ROTATION_FACE1 ;
                
                
                
                
                
            end
            
            INFOLINES.END_NODES{imat1D} =  END_NODES_GLO_NEW ;
            INFOLINES.END_ELEMENTS{imat1D} =  END_ELEMENTS_GLO_NEW ;
            
    end
end
