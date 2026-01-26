function  [DOFr,dR,NODES_LINES_ROTATIONS] =  DirichletBCglobalRVE(NODES_LINES,ndim,DISPLOC,MESH2D,DATA_REFMESH,V,DATAIN)  ;

if nargin == 0
    load('tmp.mat')
    DISPLOC = {1;0;0;0;0;0} ; 
    
    
end

DISPLOC = cell2mat(DISPLOC) ;  % Global displacement 
DISP_LOC_ALL = []; 

   ROT_RELATIVE = DATA_REFMESH.RotationMatrixFace  ;  % Rotation matrix of each of the faces of the  RVE (relative to the domain )
   
uPRESCRIBED =[] ;   

NODES_LINES_ROTATIONS = cell(length(NODES_LINES),1) ; 

MESH2D = DefaultField(MESH2D,'rotDOM',[]) ; 

for inode = 1:length(NODES_LINES)
    
    
    % Determination rotation matrix of each INTERFACE (global)
    [eDOM,iNODEloc] =  find(NODES_LINES(inode) ==  MESH2D.CN) ;
    
    eDOM = eDOM(1) ; % Number of domains containing the node
    iface = iNODEloc(1) ;  % FACE 
    if  isempty(MESH2D.rotDOM)
        Rot_GLO_INTFloc = [] ; 
    else
    Rot_GLO_DOMloc = MESH2D.rotDOM{eDOM} ;  % Rotation matrix for the domain
    Rot_DOMloc_INTFloc = ROT_RELATIVE{iface}  ;
    if ~isempty(Rot_DOMloc_INTFloc)
        Rot_GLO_INTFloc = Rot_GLO_DOMloc*Rot_DOMloc_INTFloc ;
    else
        Rot_GLO_INTFloc = Rot_GLO_DOMloc ;
    end
    end
    
    NODES_LINES_ROTATIONS{inode} =  Rot_GLO_INTFloc ; 
    
    % Then    Rot_GLO_INTFloc*V*a_GLO = V*a_LOC --> We have to solve for
    % a_LOC --< Least-squares (using geometric  mass matrix )
    Vloc = V{iface} ;   % Basis matrix (interface modes )
%     if DATAIN.ndimSP ==3 
%         nRB = 6 ;
%     else
%         nRB = 3; 
%     end
    % Nodal global displacement  ---in global coordinates
    Vloc = Vloc(:,1:length(DISPLOC)) ; 
    u = Vloc*DISPLOC ; 
    % Transformation to local coordinates
    u = RotateMatrix(Rot_GLO_INTFloc',u) ; 
   
    % Mass matrix 
    M1d = DATA_REFMESH.GeometricMassMatrixInterface{iface} ;
    ndimLOC = DATAIN.ndimSP ;
    Mbar  = sparse(size(M1d,1)*ndimLOC,size(M1d,2)*ndimLOC) ;
    for idim = 1:ndimLOC
        Mbar(idim:ndimLOC:end,idim:ndimLOC:end) = M1d ;
    end
    % Least-squares problem
    uBARloc = (Vloc'*Mbar*Vloc)\(Vloc'*Mbar*u) ;   
    
    
    uPRESCRIBED = [uPRESCRIBED;uBARloc ] ; 
    
    
end







% Old method (before 16-July-2019)
DOFS = small2large(NODES_LINES,ndim) ; % Numbering of DOFs assuming ndim = ndimMAX
nnodes = length(NODES_LINES) ;
% Known DOFs
dR = [] ;
DOFr = [] ;
ndimLOC = length(DISPLOC) ;
PRESCRIBED_ALL = 0 ;
for idim = 1:ndim
    if idim <=ndimLOC
    
            DOFr = [DOFr; (idim:ndim:length(DOFS))'] ;
            dR = [dR ; uPRESCRIBED(idim:ndimLOC:end)] ;
            PRESCRIBED_ALL = 1;
            % warning('Set it to 1')
     
    elseif idim > ndimLOC
        if PRESCRIBED_ALL == 1
            DOFr = [DOFr; (idim:ndim:length(DOFS))'] ;
            dR = [dR ; 0*ones(nnodes,1)] ;
        end
    end
end
DOFr = DOFS(DOFr) ;