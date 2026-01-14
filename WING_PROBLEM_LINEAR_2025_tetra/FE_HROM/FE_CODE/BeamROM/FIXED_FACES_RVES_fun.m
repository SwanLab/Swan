function [Gb,dR,DOFr,DOFm,NODES_ENTITIES,ANG_ROTATION_TOTAL] = ...
    FIXED_FACES_RVES_fun(a,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0
    load('tmp.mat')
end


DOFm = [] ;
DOFr = [] ;
AREA = [] ;
Gb = [] ;
dR = [] ;



% ----------------
% FACES 1 and 3 **
% ----------------

ifaces = [1,3] ;
idomxs ={1,size(DOMAINVAR.NODES_faces,1) } ;
idomys ={1:size(DOMAINVAR.NODES_faces,2),1:size(DOMAINVAR.NODES_faces,2)   } ;
NODES_LINES = {} ; 
ANG_ROTATION_TOTAL = cell(1,4) ; 
DATA = DefaultField(DATA,'angROTATION_FACE',[]) ;
ROT_MATRIX_LOCAL = cell(1,2) ;
if ~isempty(DATA.angROTATION_FACE)
    % Rotation around the y-axis
    % The normal to face 3 subtends an angle ndomx*DATA.angROTATION_FACE
    % with the global-x axis
    ndomX = size(DOMAINVAR.ListElements,1) ;
    total_angle = ndomX*DATA.angROTATION_FACE ;
    % Therefore, the rotation matrix is given by
    ROTMATRIX_f3 = [cos(total_angle),0,sin(total_angle)
        0                1    0
        -sin(total_angle)  0   cos(total_angle)] ;
    ROT_MATRIX_LOCAL = {[],ROTMATRIX_f3} ;
end

[DOFr12,dR12,nodesFA,nodesFB] =  LocPrescribedDisplac(ifaces,idomxs,idomys,DOMAINVAR,COOR,CONNECTb,TypeElementB,a,ROT_MATRIX_LOCAL)  ;
NODES_ENTITIES{1} = nodesFA; 
NODES_ENTITIES{3} = nodesFB;
ANG_ROTATION_TOTAL{1} = ROT_MATRIX_LOCAL{1} ; 
ANG_ROTATION_TOTAL{3} = ROT_MATRIX_LOCAL{2} ; 

% ----------------
% FACES 2 and 4 **
% ----------------
ROT_MATRIX_LOCAL = {[],[]} ;
ifaces = [2,4] ;
idomxs ={1:size(DOMAINVAR.NODES_faces,1), 1:size(DOMAINVAR.NODES_faces,1)} ;
idomys ={1,size(DOMAINVAR.NODES_faces,2)   } ;
[DOFr34,dR34,nodesFA,nodesFB] =  LocPrescribedDisplac(ifaces,idomxs,idomys,DOMAINVAR,COOR,CONNECTb,TypeElementB,a,ROT_MATRIX_LOCAL)  ;


NODES_ENTITIES{2} = nodesFA; 
NODES_ENTITIES{4} = nodesFB;

DOFr = [DOFr12; DOFr34] ;
dR = [dR12; dR34] ;

end

function  [DOFr,dR,nodesfA,nodesfB] =  LocPrescribedDisplac(ifaces,idomxs,idomys,DOMAINVAR,COOR,CONNECTb,TypeElementB,a,ROT_MATRIX_LOCAL)


ndim = size(COOR,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FACE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DOFr = [] ;
dR = [] ;
iface=ifaces(1) ; % Index face --
idomx =idomxs{1} ;  % Plane x= xMIN
idomy =  idomys{1} ;
nodesfLOC = (cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;


% The set of nodes may be repeated. We remove in a consistent fashion the
% repeated nodes
% --------------------------------
N = nodesfLOC(:) ;
[~,IDX ]= unique(N) ;
nodesfA = N(sort(IDX)) ;


%nodesfA = nodesfA(:) ;
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; %
% Local connectivities
% ----------------------
LocalConnect =  (CONNECTb(idomx,idomy,iface)) ;
if size(LocalConnect,1) ==1
    LocalConnect = LocalConnect' ;
end
LocalConnect = cell2mat(LocalConnect) ;
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,LocalConnect,TypeElementB) ;

COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid

if ~isempty(ROT_MATRIX_LOCAL{1})
    %   Rotation of global coordinates to local coordinates
    COORrelA = (ROT_MATRIX_LOCAL{1}'*COORrelA')' ;
end

R = ConstructBasisRigidBody(COORrelA) ;  % Rigid body modes on LOCAL COORDINATES

%%%
 

[UnrestrainedDOFS,a] = UnrestrainedDOFSlocal(a,iface,DOFA,ndim) ;



if  ~isempty(a{iface})   
    dRLOC =     R*a{iface}(:) ;   % Displacement in local coordinates    
    if ~isempty(ROT_MATRIX_LOCAL{1})
        % Rotation to global coordinates
        dRLOC = reshape(dRLOC,ndim,[]) ;
        dRLOC = ROT_MATRIX_LOCAL{1}*dRLOC ;
    end    
    
    
       dRLOC(UnrestrainedDOFS) = [] ; 
    DOFA(UnrestrainedDOFS) = [] ; 
    
    
    dR = [dR ; dRLOC(:)] ;
    DOFr = [DOFr; DOFA] ;
end


% FACE 2
% -------
iface=ifaces(2) ; % Index face --
idomx =idomxs{2} ;  % Plane x= xMIN
idomy =  idomys{2} ;
nodesfLOC = (cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
N = nodesfLOC(:) ;
[~,IDX ]= unique(N) ;
nodesfB = N(sort(IDX)) ;
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B

% if  ~isempty(a{iface})
% dR = [dR ; R*a{iface}(:)] ;
% DOFr = [DOFr; DOFB] ;
% end
%%%
[UnrestrainedDOFS,a] = UnrestrainedDOFSlocal(a,iface,DOFB,ndim) ;   % 9-July-2019

if  ~isempty(a{iface})
    
    dRLOC =     R*a{iface}(:) ;  % Local coordinates
    
    if ~isempty(ROT_MATRIX_LOCAL{2})
        dRLOC = reshape(dRLOC,ndim,[]) ;
        dRLOC = ROT_MATRIX_LOCAL{2}*dRLOC ;  % Rotation to global coordinates
        
    end
    
    dRLOC(UnrestrainedDOFS) = [] ; 
    DOFB(UnrestrainedDOFS) = [] ; 
    
    dR = [dR ; dRLOC(:)] ;
    DOFr = [DOFr; DOFB] ;
end

% Rorth = R'*M --->  R(INDX)'*N'*W*N

% Rorth = zeros(size(R')) ;
% for idim =1:3
%     INDLOC =idim:3:size(R,1) ;
%     Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ;
% end

%
% [Gb,dR,DOFr,DOFm] = BCs_BEAMS_PERIODIC(DOFA,DOFB,R,a_A,a_B,Rorth') ;

end