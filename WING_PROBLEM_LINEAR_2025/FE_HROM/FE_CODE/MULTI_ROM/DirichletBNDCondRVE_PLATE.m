function [DOFr,DOFl,dR] = DirichletBNDCondRVE_PLATE(DATAROM,MESH2D,DISP,DATAIN,ndim,DATA_REFMESH)

if nargin == 0
    load('tmp0.mat')
end

nnode = size(MESH2D.COOR,1) ; % Number of 2D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements
nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)


NODES_LINESall = MESH2D.NODES_LINESall ;  % Both midside and corner points
nlines = min(length(NODES_LINESall),length(DISP.LINE)) ;


% ----------
dR = [] ;
DOFr = [] ;

for ilines = 1:nlines
    
    [xC,yC,L_total,CNb] = CoordLineElement(NODES_LINESall,MESH2D,ilines) ;
    % Now we select corner points
    NodVert= intersect(MESH2D.NODESvert,CNb(:));
    NODES_LINES = zeros(size(NodVert)) ;
    % New numbering
    for  inode = 1:length(NODES_LINES)
        NODES_LINES(inode) =  find(NodVert(inode)==MESH2D.NODESvert) ;
    end
    % Coordinates
    COORloc = MESH2D.COORvert(NODES_LINES,:) ;
    % Coordinates relative to the centroid
    COORrel = bsxfun(@minus,COORloc',[xC,yC,0]')' ;
    % Rigid  Body modes
    R = ConstructBasisRigidBody(COORrel) ;
    
    [DOFrLOC,dRloc] = BCsSingleLine(R,NODES_LINES,ilines,DISP,ndim,COORrel) ;
    
    DOFr = [DOFr; DOFrLOC] ;
    dR = [dR; dRloc] ;
end

% Remove repeated DOFs
[DOFr,IIII] = unique(DOFr) ; 
dR = dR(IIII) ; 

 
DOFl = 1:ndim*nnode ;
DOFl(DOFr) = [] ;




end
%
function  [DOFr,dR] = BCsSingleLine(R,NODES_LINES,ilines,DISP,ndim,COORrel)

DOFS = small2large(NODES_LINES,ndim) ; % All DOFs associated to the nodes of the studied line (x 5 
%). These are possible candidates for being prescribed DOFs 
nnodes = length(NODES_LINES) ;
DISPLOC  =DISP.LINE{ilines} ; 
dR = zeros(size(DOFS)) ; 
VISITED = zeros(size(DOFS)) ; 


ndimLOC = 3 ; 
% Translations  
% ------------------------------ 
for  imode  =1:ndimLOC
    DISPLOC = DISP.LINE{ilines}(imode) ;
    idim = imode;
    if ~isempty(DISPLOC{1})
        dRloc = R(idim:ndimLOC:end,idim)*DISPLOC{1} ;  % Imposed displacements
        DOFrLOC = idim:ndim:length(DOFS) ;  % Associated DOFs
        dR(DOFrLOC) =  dR(DOFrLOC) + dRloc ;
        VISITED(DOFrLOC) =  1;
    end
end

% Rotations in the x,y,z-direction
% ---------------------------
for imode = 4:6
    DISPLOC = DISP.LINE{ilines}(imode) ;
    if ~isempty(DISPLOC{1})
        % Rotation of the corner line
        if  imode ~=6
        idim = imode ;
        DOFrLOC = idim:ndim:length(DOFS) ;
        dRloc = ones(size(DOFrLOC))*DISPLOC{1} ;
        dR(DOFrLOC) =   dRloc ;
        VISITED(DOFrLOC) =  1;
        end
        % Associated translations (we assume all are prescribed)
        idimGLO = [1:3] ;
        for idimLOC = 1:length(idimGLO)
            idim = idimGLO(idimLOC) ;
            DOFrLOC = idim:ndim:length(DOFS) ;
            dRloc = R(idim:ndimLOC:end,imode)*DISPLOC{1} ;
            dR(DOFrLOC) =  dR(DOFrLOC) + dRloc ;
            VISITED(DOFrLOC) =  1;
        end
        
    end
end

SELECTED_DOFS = find(VISITED==1) ; 

dR = dR(SELECTED_DOFS) ; 
DOFr = DOFS(SELECTED_DOFS)  ; 

 



 


% 
% % Known DOFs
% dR = [] ;
% DOFr = [] ;
% 
% for idim = 1:ndim
%     if ~isempty(DISPLOC{idim})
%         DOFr = [DOFr; (idim:ndim:length(DOFS))'] ;
%         dR = [dR ; DISPLOC{idim}*ones(nnodes,1)] ;
%     end
% end
% DOFr = DOFS(DOFr) ;

end
%
%
%
% %
% %
% % % RIGHT END
% % %----------------
% %
% %
% % NODE_2 = MESH2D.RIGHT_END_NODE ;
% % DOFS_2= small2large(NODE_2,ndim) ;
% % % Known DOFs
% % r = [] ;
% % DISPLOC  =DISP.RIGHT_END ;
% % for idim = 1:length(DISPLOC)
% %     if ~isempty(DISPLOC{idim})
% %         r = [r; idim] ;
% %         dR = [dR ; DISPLOC{idim}] ;
% %     end
% % end
% % DOFr = [DOFr; DOFS_2(r)] ;
% %
% %
% %
% % DOFl = setdiff(1:ndim*nnode,DOFr) ;
% %
% %
