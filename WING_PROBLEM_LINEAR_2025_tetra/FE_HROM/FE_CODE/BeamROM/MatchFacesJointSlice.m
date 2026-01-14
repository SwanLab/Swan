function [FACES_ORDER, ROTATION_FACE1,SUCCESSFUL_MATCH ]= ...
    MatchFacesJointSlice(CENTROIDslices,COORslices,CENTROIDjoints,COORjoint,TOL_ABSOLUTE)
%
%  Given
%  CENTROIDslices{i} --> Coordinates rotated centroid of the SLICE
%  COORslices{i}     ---> Coordinates rotated nodal faces of the SLICE
%  CENTROIDjoints{j} --> Coordinates reference centroid of the JOINT
%  COORjoint{j}    --> Coordinates reference nodal faces  of the JOINT
%
%
%  MatchFacesJointSlice returns the correspondance between rotated (SLICES) and
%  reference objects (JOINTS), as well as the corresponding rotation (with respect the centroid of FACE1)
%
%  ------
%  JAHO, 8-FEb-2018, revised 29-May-2018, 11-July-2018
%dbstop('16')
if nargin == 0
    load('tmp1.mat')
end
SUCCESSFUL_MATCH = 1;
ROTATION_FACE1 =zeros(3);
FACES_ORDER  = 1:length(COORslices) ;
TOL_LENGTHS = TOL_ABSOLUTE;
% We are going to determine which SLICE face corresponds to the reference face = 1
%---------------------------------------------------------------------
ifaceJOINT = 0;
% But what happens if face 1 has no slice ---for instance, it corresponds
% to a built-end ?
% if isempty(COORslices{ifaceJOINT})  % If the first one is empty
%     ifaceJOINT = 2 ;
% end

% We wish to determine which entry of COORslices  (GLOBAL NUMBERING) corresponds to COORjoint{ifaceJOINT}
% (LOCAL NUMBERING)
% To this end, we begin by picking up  3 non-aligned points of the first face of the
% reference geometry under study.  (of the joint )
TOL = 1e-16 ;
FOUND = 0 ;
while ifaceJOINT < length(COORjoint)  && FOUND == 0
    ifaceJOINT = ifaceJOINT+1 ;
    [Xref NODES] = licols(COORjoint{ifaceJOINT},TOL); % Select linearly independent columns
    NODES = NODES(1:3) ;
    Xref = Xref(:,1:3) ;
    ifaceSLICE = 1; % First guess is that the rotated face is the first one
    FOUND = 0 ;
    while ifaceSLICE <= length(COORslices)
        % First check:  is   the number of nodes the  same ?
        if size(COORjoint{ifaceJOINT},2) ~= size(COORslices{ifaceSLICE},2)
            ifaceSLICE = ifaceSLICE + 1;
        else
            % Second check: since nodes are already ordered in the joint so that they  match
            % those of the slices, if the matching fails, it means that this
            % is not the sought face
            % ------------------------------------------------------
            Xrot =  COORslices{ifaceSLICE}(:,NODES);  % Coordinates of the rotated configuration, slice
            [TRANSLATION,  ROTATION ]= RotationMatrixFrom3points(Xref,Xrot,TOL_LENGTHS)  ;
            if ~isempty(ROTATION)
                %   FACES_ORDER(1) = ifaceSLICE ;
                
                %  The match is almost done: same number of nodes, and
                %  correspondence between chosen nodes. Now let us check
                %  whether the remaining ends ( centroids) match each other
                % --------------------------------------------------------
                % To this end, we translate and rotate the centroids of the reference configuration
                % NOTE: Rotation has been calculated with respect  to the first
                % node of
                % set "NODES"
                ROTcentroidsREF = zeros(size(CENTROIDjoints)) ;
                %%%%%%%%%%%%%%%%%%%
                for inodeL = 1:size(CENTROIDjoints,2)
                    ROTcentroidsREF(:,inodeL) =  Xref(:,1) + TRANSLATION + ROTATION*(CENTROIDjoints(:,inodeL)-Xref(:,1)) ;
                end
                %% Now we have to compare ROTcentroidsREF (centroid joint ) with CENTROIDslices (centroid slices)
                % Let us exploit the capabilities of knnsearch
                %    [IDX,DIST] = knnsearch(ROTcentroidsREF',CENTROIDslices') ;
                [IDX,DIST] = knnsearch(CENTROIDslices',ROTcentroidsREF') ;
                if any(DIST > TOL_ABSOLUTE)
                    % This is not the sought-after face
                    ifaceSLICE = ifaceSLICE + 1 ;
                else
                    % We found the correspondence !!!
                    FACES_ORDER= IDX ;
                    FOUND = 1;
                    break
                end
                
            else
                % Rotation is not possible because the relative distance
                % bewteen the chosen points changes from one configuration to
                % the the other
                ifaceSLICE = ifaceSLICE + 1;
            end
        end
    end
    
    %  ifaceSLICE = 1;
end

if FOUND == 0
    warning('This joint does not fit in the prescribed location !!!! ')
    SUCCESSFUL_MATCH = 0;
else
    % Rotation and translation with respect the centroid of  each face 1
    % -------------------------------------------------------------------
    %  for ifaceJOINT =1:length(COORjoint)
    % ifaceJOINT  = 1 ;
    %    [Xref NODES] = licols(COORjoint{ifaceJOINT}); %
    % REFERENCE
    XREF1 =   CENTROIDjoints(:,ifaceJOINT) ;
    NODE2 = NODES(1) ;
    XREF2 = COORjoint{ifaceJOINT}(:,NODE2) ;
    NODE3 = NODES(2) ;
    XREF3 = COORjoint{ifaceJOINT}(:,NODE3) ;
    XREF = [XREF1,XREF2,XREF3] ;
    if rank(XREF) ~=3
        NODE3 = NODES(3) ;
        XREF3 = COORjoint{ifaceJOINT}(:,NODE3) ;
        XREF = [XREF1,XREF2,XREF3] ;
    end
    % ROTATED
    ifaceSLICE = FACES_ORDER(ifaceJOINT) ;
    XR1 =   CENTROIDslices(:,ifaceSLICE) ;
    XR2 = COORslices{ifaceSLICE}(:,NODE2) ;
    XR3 = COORslices{ifaceSLICE}(:,NODE3) ;
    XROT = [XR1,XR2,XR3] ;
    %
    [TRANSLATION,  ROTATION_FACE1 ]= RotationMatrixFrom3points(XREF,XROT,TOL_LENGTHS)  ;
    %   end
    
end
