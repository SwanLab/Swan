function CONNECTslices = ConnectingSlices_JOINT(END_NODES,END_ELEMENTS,MESH1D,MESH3D)
if nargin == 0
    load('tmp1.mat')
end
% See BeamROM.pdf
% Coordinates of the 3D nodes of the slices connected with joints whose end
% nodes and end elements are END_NODES and END_ELEMENTS
% OUTPUT
%CONNECTslices.COOR = NODES_FACES_SLICES_COOR ;
%CONNECTslices.NODES = NODES_FACES_SLICES ;
%CONNECTslices.TRANSLATION = TranslationVectorGLO ;  --> This is just for
%GID batch file
% JAHO, 6-Feb-2018
% -------------------------------
ndim = size(MESH1D.COOR,2) ;
for inodeLOC =  1:length(END_NODES)
    inode = END_NODES(inodeLOC) ; % Loop over end nodes of the joint
    [ielem,dummy] = find(MESH1D.CN == inode) ;
    ielem = setdiff(ielem,END_ELEMENTS(inodeLOC)) ; % Number of element corresponding to end node "inode"
    
    if isempty(ielem)
           NODES_FACES_SLICES_COOR{inodeLOC} = [] ;
        NODES_FACES_SLICES{inodeLOC} = [] ;
        ROTATION_FACES_SLICES{inodeLOC} = [] ;
       
        
        COORtransfGLO{inodeLOC} = [] ;
        TranslationVectorGLO{inodeLOC} = [] ;
        indexSLICE_GLO(inodeLOC) = 0 ;
        
    else
        % Type of beam associated to elem "ielem"
        itypeBEAM = MESH1D.MaterialType(ielem) ;
        % Subtype of beam
        iBEAM = MESH1D.INFOLINES.SUBTYPESentity(ielem) ;
        % Nodes comprising the beam
        [NODESbeam] = MESH1D.INFOLINES.NODES{itypeBEAM}{iBEAM} ;
        % Is inode the final or starting point of the beam ?
        PROPS =  MESH1D.PROP(itypeBEAM) ;
        PROPS = DefaultField(PROPS,'xLOCAL',1) ; % = 1 ;
        SIGNO = PROPS.xLOCAL  ;
        if isempty(SIGNO) ; SIGNO = 1 ; end
        if SIGNO == -1
            NODESbeam = NODESbeam(end:-1:1) ;
        end
        xREF = MESH1D.COOR(inode,:)  ;
        if inode == NODESbeam(1)
            ISstart =1;
            % Coordinates of the initial point (face 1), centroid
            xINI = MESH1D.COOR(inode,:) ;
        else
            ISstart =0;
            inodeINI = setdiff(MESH1D.CN(ielem,:),inode) ;
            xINI = MESH1D.COOR(inodeINI,:) ;
        end
        
        
        % Corresponding slice
        indexSLICE = MESH1D.PROP(itypeBEAM).INDEX_ENTITY ;
        % For mixed beams, we take the first slice
        indexSLICE = indexSLICE(1) ;
        % Rotation matrix
        ifin = ielem*ndim ; iini = ielem*ndim-ndim+1 ;
        R = MESH1D.ROTATIONS(:,iini:ifin) ;
        % Transformation matrices
        a0 = MESH1D.TRANSFM.a0(:,ielem) ;
        A = MESH1D.TRANSFM.A(:,iini:ifin)  ;
        D = MESH1D.TRANSFM.D(:,iini:ifin)  ;
        % NODES 3D mesh in contact with the joint
        iface =1  ; % B
        CENTROIDopp =  MESH3D.SLICES(indexSLICE).DATA3D.CENTRf2  ;
        if ISstart == 0
            iface = 2 ;
            CENTROIDopp =  MESH3D.SLICES(indexSLICE).DATA3D.CENTRf1 ;
        end
        NODESface  = MESH3D.SLICES(indexSLICE).DATA3D.NODES_FACES{iface} ;
        % Transformed coordinates
        COORrel  = zeros(size(MESH3D.SLICES(indexSLICE).DATA3D.COOR')) ;
        for idim=1:ndim
            COORrel(idim,:) = MESH3D.SLICES(indexSLICE).DATA3D.COOR(:,idim) - MESH3D.SLICES(indexSLICE).DATA3D.CENTRf1(idim) ;
        end
        % Transformation (scaling, varying cross-section)
        % x = a0 + A*X +X_1*D*X
        COORtransf = zeros(size(COORrel)) ;
        for idim = 1:ndim
            COORtransf(idim,:) = a0(idim) + A(idim,:)*COORrel + COORrel(1,:).*(D(idim,:)*COORrel) ;
        end
        % Rotated coordinates
        COORrel = R*COORtransf ;
        % Translation
        for idim = 1:ndim
            COORrel(idim,:) = COORrel(idim,:) + xINI(idim) ;
        end
        % -----------------------------------------------------------------------
        TranslationVector =xREF'- CENTROIDopp  ;
        % -----------------------------------------------------------------------
        % Coordinate face
        % ---------------
        COORtransf = COORrel(:,NODESface) ;
        
        NODES_FACES_SLICES_COOR{inodeLOC} = COORtransf ;
        NODES_FACES_SLICES{inodeLOC} = NODESface ;
        ROTATION_FACES_SLICES{inodeLOC} = R ;
        CENTROID_transf = sum(COORtransf,2)/size(COORtransf,2) ; %
        
        COORtransfGLO{inodeLOC} = COORtransf ;
        TranslationVectorGLO{inodeLOC} = TranslationVector ;
        indexSLICE_GLO(inodeLOC) = indexSLICE ;
    end
    
    
end


CONNECTslices.COOR = NODES_FACES_SLICES_COOR ;
CONNECTslices.NODES = NODES_FACES_SLICES ;
CONNECTslices.TRANSLATION = TranslationVectorGLO ;
CONNECTslices.ROTATION = ROTATION_FACES_SLICES ;
CONNECTslices.indexSLICE_GLO = indexSLICE_GLO ;


