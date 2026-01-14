function [Gb,dR,DOFr,DOFm] = PRESCRIBED_DISP_PLATES_RVE_fun(PRESCRIBED_VALUES,NAMEPROJECT,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)
% Imposition of BCs on plates,  JAHO, 1-Sept-2018
% See also LocalNamesPRESCRIBED_plates.m 
if nargin == 0
    load('tmp1.mat')
end

ndim = 3;
DOFm = [] ;
DOFr = [] ;
Gb = [] ;
dR = [] ;

[x,y,z,COOR_FACE,nodesf] =  NodesFacesCorners(DOMAINVAR,COOR)  ;

% Determination of corner nodes  and side nodes
% ----------------------------------------------
DATAINP.NODES_FACES  = {nodesf.A,nodesf.B,nodesf.C,nodesf.D} ; 
[NODES_CORNERS,NODES_SIDES] = CornerSideNodes(DATAINP) ; 

% Boundary nodes (corner1,...corner4,side1,...side4)
% --------------------------------------------------
BNODES = [cell2mat(NODES_CORNERS');cell2mat(NODES_SIDES')] ;  
% Q: Matrix of global modes of the domain (NDOFboundaryx20). DOFboundary = [corners, sides]  ; 
% % Modes corresponding to assuming tri-linear state of displacement 
% -----------------------------------------------------------------------------------------------
Q = GlobalModesPlate(COOR,BNODES) ; 
% Matrix relating global modes amplitude with amplitudes of corner modes 
S = NodalValuesCornersToGlobalModes(NODES_CORNERS,x,y,z,COOR,Q) ; 
% Amplitude modes for tests 
Ly = y.max - y.min ;  
Lx = x.max - x.min ; 
% -------------------
% Set of 14 modes (linear dependence on x,y,z)
% --------------------
[aTESTS,NAMES ]= DesignTestPlates(PRESCRIBED_VALUES.DISPLACEMENT,PRESCRIBED_VALUES.ROTATION,Ly,Lx) ; 
% --------------------------------------------------------------------------
% Check which kind of test is going to be conducted
% ---------------------------------------------------------------------------

USE_QUADRATIC_BENDING_MODES = 1 ;  % Rather than linear modes 

if USE_QUADRATIC_BENDING_MODES == 1 
    DOFr = small2large(BNODES,ndim) ; 
    k = 0.9 ; 
    COORB = COOR(BNODES,:) ;
    ALPHA = PRESCRIBED_VALUES.ROTATION ; 
    % Centroid 
    xC(1) = 0.5*(x.max + x.min) ; 
    xC(2) = 0.5*(y.max + y.min) ; 
    xC(3) = 0.5*(z.max + z.min) ; 
    COORB = bsxfun(@minus,COORB',xC')';
     dR = zeros(length(DOFr),1) ;
    switch NAMEPROJECT
        case {'BendingXZ'}
              
            idim = 1; 
            idir = 1; 
            L = Lx ; 
            dR(idim:ndim:end) = 2*ALPHA/L*COORB(:,idir).*COORB(:,3) ; ;
            %%%
            idim = 3; 
            dR(idim:ndim:end) =   ALPHA*L/4 - ALPHA/L*COORB(:,idir).^2 ;
            
        case {'BendingYZ'}
            
             idim = 2; 
            idir = 2; 
            L = Ly ; 
            dR(idim:ndim:end) = 2*ALPHA/L*COORB(:,idir).*COORB(:,3) ; ;
            %%%
            idim = 3; 
            dR(idim:ndim:end) =   ALPHA*L/4 - ALPHA/L*COORB(:,idir).^2 ;
            
                
        case {'StrainXZ'}
           
            idim = 1; 
            idir = 1; 
            L = Lx ; 
            dR(idim:ndim:end) = (ALPHA*COORB(:,3)/(2*L)).*(L + 2*COORB(:,idir)) ; 
            %%%
            idim = 3; 
            dR(idim:ndim:end) =  -ALPHA*(L + 2*COORB(:,idir) ).*(2*COORB(:,idir) -L + 2*L*k)/(8*L) ;
                    
        case {'StrainYZ'}
             idim = 2; 
            idir = 2; 
            L = Ly ; 
            dR(idim:ndim:end) = (ALPHA*COORB(:,3)/(2*L)).*(L + 2*COORB(:,idir)) ; 
            %%%
            idim = 3; 
            dR(idim:ndim:end) =  -ALPHA*(L + 2*COORB(:,idir) ).*(2*COORB(:,idir) -L + 2*L*k)/(8*L) ;
            
            
        otherwise
             [dR,DOFr] = LocalDisplacements(NAMES,NAMEPROJECT,S,aTESTS,Q,BNODES,ndim) ; 
    end
    
else
     [dR,DOFr] = LocalDisplacements(NAMES,NAMEPROJECT,S,aTESTS,Q,BNODES,ndim)  ; 
    
end





% ------------------------------
% 
% % Designing the set of experiments 
% % ---------------------------------
% % Face 3 
% % ------
%  
%  
% 
% 
% 
% 
% 
% % 
% % % Boundary DOFs --->
% % NODES = [nodesf.A,nodesf]
% %  
% % DOFr = [] ;
% % dR = [] ;
% % iface=ifaces(1) ; % Index face --
% % idomx =idomxs{1} ;  % Plane x= xMIN
% % idomy =  idomys{1} ;
% % % nodesfA = nodesf.(cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
% % % nodesfA = nodesfA(:) ;
% % nodesfA = nodesf.(FACESNAMES{iface}) ; 
% % DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
% % COOR_FACE = COOR(nodesfA,:) ; %
% % 
% % 
% 
% % 
% % 
% % 
% % %%%
% % dR = [dR ; R*a{iface}'] ;
% % DOFr = [DOFr; DOFA] ;
% 
% 
% 
% %%%
% if any((a{3}))
%     [DOFr,dR] = PrescribedDisplacement_FACE_C(a,DOMAINVAR,x,y,z,COOR,CONNECTb,TypeElementB,COOR_FACE,nodesf,ndim) ;
% elseif any(a{4})
%     [DOFr,dR] = PrescribedDisplacement_FACE_D(a,DOMAINVAR,x,y,z,COOR,CONNECTb,TypeElementB,COOR_FACE,nodesf,ndim) ;
% elseif  any(a{2})
%     [DOFr,dR] = PrescribedDisplacement_FACE_D(a,DOMAINVAR,x,y,z,COOR,CONNECTb,TypeElementB,COOR_FACE,nodesf,ndim) ;
% else  
%     error('Option not implemented')
% end
% 
% 
end

    function [dR,DOFr] = LocalDisplacements(NAMES,NAMEPROJECT,S,aTESTS,Q,BNODES,ndim)
        
        
        
        icheck = 1; EXIT = 0 ;
        while icheck <= length(NAMES) & EXIT == 0
            ok = strcmp(NAMES{icheck},NAMEPROJECT) ;
            if ok ==1
                EXIT = 1;
            else
                icheck = icheck+1 ;
            end
        end
        
        if EXIT == 0
            error('Project not found')
        else
            iproj = icheck  ;
        end
        
        aCORNER = aTESTS(:,iproj) ;
        % Therefore, the set of amplitudes for the nodal modes of the boundary is:
        qNODES = S*aCORNER;
        
        % Therefore
        dR = Q*qNODES ;
        DOFr = small2large(BNODES,ndim) ;
        
        
    end



% 
% function  [DOFr,dR,dP] =  LocPrescribedDisplacCOOR(ifaces,idomxs,idomys,DOMAINVAR,COOR,CONNECTb,TypeElementB,a,P,nodesf)
% 
% 
% ndim = 3;
% 
% FACESNAMES =fieldnames(nodesf) ; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% FACE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOFr = [] ;
% dR = [] ;
% iface=ifaces(1) ; % Index face --
% idomx =idomxs{1} ;  % Plane x= xMIN
% idomy =  idomys{1} ;
% % nodesfA = nodesf.(cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
% % nodesfA = nodesfA(:) ;
% nodesfA = nodesf.(FACESNAMES{iface}) ; 
% DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
% COOR_FACE = COOR(nodesfA,:) ; %
% 
% 
% % Local connectivities
% % ----------------------
% LocalConnect =  (CONNECTb(idomx,idomy,iface)) ;
% if size(LocalConnect,1) ==1
%     LocalConnect = LocalConnect' ;
% end
% LocalConnect = cell2mat(LocalConnect) ;
% [CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,LocalConnect,TypeElementB) ;
% 
% COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
% R = ConstructBasisRigidBody(COORrelA) ;
% 
% 
% 
% %%%
% dR = [dR ; R*a{iface}'] ;
% DOFr = [DOFr; DOFA] ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Displacement corner points
% % ---------------------------
% X = P.AB ;
% for i = 1:2
%     Xrel = X{i} -  CentroidFA ;
%     Rloc = ConstructBasisRigidBody(Xrel) ;
%     dP.AB{i} = Rloc*a{iface}' ;
% end
% X = P.AD ;
% for i = 1:2
%     Xrel = X{i} -  CentroidFA ;
%     Rloc = ConstructBasisRigidBody(Xrel) ;
%     dP.AD{i} = Rloc*a{iface}' ;
% end
% 
% 
% 
% 
% 
% % FACE 2
% % -------
% iface=ifaces(2) ; % Index face --
% idomx =idomxs{2} ;  % Plane x= xMIN
% idomy =  idomys{2} ;
% % nodesfB = (cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
% % nodesfB = nodesfB(:) ;
% nodesfB = nodesf.(FACESNAMES{iface}) ; 
% DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B
% 
% COOR_FACE = COOR(nodesfB,:) ; %
% COOR_FACE_B = COOR_FACE ;
% 
% dR = [dR ; R*a{iface}'] ;
% DOFr = [DOFr; DOFB] ;
% 
% 
% % ---------------------------
% X = P.CB ; %  
% for i = 1:2
%     Xrel = X{i} -  CentroidFA ;
%     Rloc = ConstructBasisRigidBody(Xrel) ;
%     dP.CB{i} = Rloc*a{iface}' ;
% end
% X = P.CD ;
% for i = 1:2
%     Xrel = X{i} -  CentroidFA ;
%     Rloc = ConstructBasisRigidBody(Xrel) ;
%     dP.CD{i} = Rloc*a{iface}' ;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% % Rorth = R'*M --->  R(INDX)'*N'*W*N
% 
% % Rorth = zeros(size(R')) ;
% % for idim =1:3
% %     INDLOC =idim:3:size(R,1) ;
% %     Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ;
% % end
% 
% %
% % [Gb,dR,DOFr,DOFm] = BCs_BEAMS_PERIODIC(DOFA,DOFB,R,a_A,a_B,Rorth') ;
% 
% end
% 


function [x,y,z,COOR_FACE,nodesf] =  NodesFacesCorners(DOMAINVAR,COOR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORNER POINTS
%---------------------
% FACE A (face 1)
% ------
iface  = 1 ;
idomx =1 ;
idomy =1:size(DOMAINVAR.NODES_faces,2);
nodesfLOC = (cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;

% The set of nodes may be repeated. We remove in a consistent fashion the
% repeated nodes
% --------------------------------
N = nodesfLOC(:) ; 
[~,IDX ]= unique(N) ; 
nodesf.A = N(sort(IDX)) ; 

COOR_FACE.A = COOR(nodesf.A,:) ; %
x.min = COOR_FACE.A(1,1) ;
z.max = max(COOR_FACE.A(:,3)) ;
z.min = min(COOR_FACE.A(:,3)) ;


% FACE B (face 2)
% ------
iface  = 2 ;
idomy = 1;
idomx =1:size(DOMAINVAR.NODES_faces,1);
nodesfLOC = (cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
N = nodesfLOC(:) ; 
[~,IDX ]= unique(N) ; 
nodesf.B = N(sort(IDX)) ; 
COOR_FACE.B = COOR(nodesf.B,:) ; %
y.min = COOR_FACE.B(1,2) ;


% FACE C (face 3)
% ------
iface  = 3 ;
idomx = size(DOMAINVAR.NODES_faces,1);
idomy =1:size(DOMAINVAR.NODES_faces,2);
nodesfLOC = (cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
N = nodesfLOC(:) ; 
[~,IDX ]= unique(N) ; 
nodesf.C = N(sort(IDX)) ; 

COOR_FACE.C = COOR(nodesf.C,:) ; %
x.max = COOR_FACE.C(1,1) ;




% FACE D  (4)
% ------
iface  = 4;
idomy = size(DOMAINVAR.NODES_faces,2);
idomx =1:size(DOMAINVAR.NODES_faces,1);
nodesfLOC = (cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
N = nodesfLOC(:) ; 
[~,IDX ]= unique(N) ; 
nodesf.D = N(sort(IDX)) ; 
COOR_FACE.D = COOR(nodesf.D,:) ; %
y.max = COOR_FACE.D(1,2) ;

end

% 
% function   [DOFr,dR] = PrescribedDisplacement_FACE_C(a,DOMAINVAR,x,y,z,COOR,CONNECTb,TypeElementB,COOR_FACE,nodesf,ndim)
% 
% % Case in which the prescribed displacement is on face C
% if any(a{2}) || any(a{4})
%     error('Displacements should be prescribed on one single face')
% end
% % ----------------
% % FACES 1 and 3 ** (faces A and C)
% % ----------------
% ifaces = [1,3] ;
% idomxs ={1,size(DOMAINVAR.NODES_faces,1) } ;
% idomys ={1:size(DOMAINVAR.NODES_faces,2),1:size(DOMAINVAR.NODES_faces,2)   } ;
% % Corner points face A
% P.AB{1} = [x.min,y.min,z.min];
% P.AB{2} = [x.min,y.min,z.max];
% P.AD{1} = [x.min,y.max,z.min];
% P.AD{2} = [x.min,y.max,z.max];
% % Corner points face C
% P.CB{1} = [x.max,y.min,z.min];
% P.CB{2} = [x.max,y.min,z.max];
% P.CD{1} = [x.max,y.max,z.min];
% P.CD{2} = [x.max,y.max,z.max];
% 
% [DOFr13,dR13,dP] =  LocPrescribedDisplacCOOR(ifaces,idomxs,idomys,DOMAINVAR,COOR,CONNECTb,TypeElementB,a,P,nodesf)  ;
% 
% %%%% Displacements faces 2 and 4 (B and D)
% % ------------------------------
% % Transformation matrix  d = J X  --> J = d*X^{-1}
% % FACE 2 (B)
% dMAT = [dP.AB{1},dP.CB{1},dP.CB{2}] ;
% xMAT  =[P.AB{1}',P.CB{1}',P.CB{2}'] ;
% J = dMAT*inv(xMAT) ;  % Transformation matrix
% % Displacement face B
% uB = J*COOR_FACE.B' ;
% DOFr2 = small2large(nodesf.B,ndim) ;
% dR2 = uB(:) ;
% 
% % FACE 4 (D)
% dMAT = [dP.AD{1},dP.CD{1},dP.CD{2}] ;
% xMAT  =[P.AD{1}',P.CD{1}',P.CD{2}'] ;
% J = dMAT*inv(xMAT) ;  % Transformation matrix
% % Displacement face D
% uD = J*COOR_FACE.D' ;
% DOFr4 = small2large(nodesf.D,ndim) ;
% dR4 = uD(:) ;
% 
% DOFr = [DOFr13;DOFr2;DOFr4] ;
% dR = [dR13; dR2; dR4] ;
% 
% 
% end
% 
% 
% 
% function   [DOFr,dR] = PrescribedDisplacement_FACE_D(a,DOMAINVAR,x,y,z,COOR,CONNECTb,TypeElementB,COOR_FACE,nodesf,ndim)
% 
% % Case in which the prescribed displacement is on face D
% if any(a{1}) || any(a{3})
%     error('Displacements should be prescribed on one single face')
% end
% % ----------------
% % FACES 2 and 4 ** (faces B and D)
% % ----------------
% ifaces = [2,4] ;
% idomys ={1,size(DOMAINVAR.NODES_faces,2) } ;
% idomxs ={1:size(DOMAINVAR.NODES_faces,1),1:size(DOMAINVAR.NODES_faces,1)   } ;
% 
% P.BA{1} = [x.min,y.min,z.min];
% P.BA{2} = [x.min,y.min,z.max];
% P.BC{1} = [x.max,y.min,z.min];
% P.BC{2} = [x.max,y.min,z.max];
% 
% P.DA{1} = [x.min,y.max,z.min];
% P.DA{2} = [x.min,y.max,z.max];
% P.DC{1} = [x.max,y.max,z.min];
% P.DC{2} = [x.max,y.max,z.max];
% 
% [DOFr24,dR24,dP] =  LocPrescribedDisplacCOOR_D(ifaces,idomxs,idomys,DOMAINVAR,COOR,CONNECTb,TypeElementB,a,P,nodesf)  ;
% 
% %%%% Displacements faces 3 and 4
% % ------------------------------
% % Transformation matrix  d = J X  --> J = d*X^{-1}
% % FACE 1 (A)
% dMAT = [dP.BA{1},dP.DA{1},dP.DA{2}] ;
% xMAT  =[P.BA{1}',P.DA{1}',P.DA{2}'] ;
% J = dMAT*inv(xMAT) ;  % Transformation matrix
% % Displacement face C
% uA = J*COOR_FACE.A' ;
% DOFr1 = small2large(nodesf.A,ndim) ;
% dR1 = uA(:) ;
% 
% % FACE 3 (C)
% dMAT = [dP.BC{1},dP.DC{1},dP.DC{2}] ;
% xMAT  =[P.BC{1}',P.DC{1}',P.DC{2}'] ;
% J = dMAT*inv(xMAT) ;  % Transformation matrix
% % Displacement face C
% uC = J*COOR_FACE.C' ;
% DOFr3 = small2large(nodesf.C,ndim) ;
% dR3 = uC(:) ;
% 
% DOFr = [DOFr1 ; DOFr3 ;DOFr24] ;
% dR = [dR1;dR3; dR24] ;
% 
% 
% end
% 
% function  [DOFr,dR,dP] =  LocPrescribedDisplacCOOR_D(ifaces,idomxs,idomys,DOMAINVAR,COOR,CONNECTb,TypeElementB,a,P,nodesf)
% 
% 
% ndim = 3;
% 
% % 
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%% FACE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % DOFr = [] ;
% % dR = [] ;
% % iface=ifaces(1) ; % Index face --
% % idomx =idomxs{1} ;  % Plane x= xMIN
% % idomy =  idomys{1} ;
% % % nodesfA = nodesf.(cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
% % % nodesfA = nodesfA(:) ;
% % nodesfA = nodesf.(FACESNAMES{iface}) ; 
% % DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
% % COOR_FACE = COOR(nodesfA,:) ; %
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% FACE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FACESNAMES =fieldnames(nodesf) ; 
% DOFr = [] ;
% dR = [] ;
% iface=ifaces(1) ; % Index face --
% idomx =idomxs{1} ;  % Plane x= xMIN
% idomy =  idomys{1} ;
% %nodesfA = (cell2mat(DOMAINVAR.NODES_faces(idomx,idomy,iface))) ;
% %nodesfA = nodesfA(:) ;
%  nodesfA = nodesf.(FACESNAMES{iface}) ; 
% DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
% COOR_FACE = COOR(nodesfA,:) ; %
% 
% 
% % Local connectivities
% % ----------------------
% LocalConnect =  (CONNECTb(idomx,idomy,iface)) ;
% if size(LocalConnect,1) ==1
%     LocalConnect = LocalConnect' ;
% end
% LocalConnect = cell2mat(LocalConnect) ;
% [CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,LocalConnect,TypeElementB) ;
% 
% COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
% R = ConstructBasisRigidBody(COORrelA) ;
% 
% 
% 
% %%%
% dR = [dR ; R*a{iface}'] ;
% DOFr = [DOFr; DOFA] ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Displacement corner points
% % ---------------------------
% X = P.BA ;
% for i = 1:2
%     Xrel = X{i} -  CentroidFA ;
%     Rloc = ConstructBasisRigidBody(Xrel) ;
%     dP.BA{i} = Rloc*a{iface}' ;
% end
% X = P.BC ;
% for i = 1:2
%     Xrel = X{i} -  CentroidFA ;
%     Rloc = ConstructBasisRigidBody(Xrel) ;
%     dP.BC{i} = Rloc*a{iface}' ;
% end
% 
% 
% 
% 
% 
% % FACE 2
% % -------
% iface=ifaces(2) ; % Index face --
% idomx =idomxs{2} ;  % Plane x= xMIN
% idomy =  idomys{2} ;
% nodesfB = nodesf.(FACESNAMES{iface}) ; 
% DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B
% 
% COOR_FACE = COOR(nodesfB,:) ; %
% COOR_FACE_B = COOR_FACE ;
% 
% dR = [dR ; R*a{iface}'] ;
% DOFr = [DOFr; DOFB] ;
% 
% 
% % ---------------------------
% X = P.DA ; % The same as P.DA (when sustracted the centroid)
% for i = 1:2
%     Xrel = X{i} -  CentroidFA ;
%     Rloc = ConstructBasisRigidBody(Xrel) ;
%     dP.DA{i} = Rloc*a{iface}' ;
% end
% X = P.DC ;
% for i = 1:2
%     Xrel = X{i} -  CentroidFA ;
%     Rloc = ConstructBasisRigidBody(Xrel) ;
%     dP.DC{i} = Rloc*a{iface}' ;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% % Rorth = R'*M --->  R(INDX)'*N'*W*N
% 
% % Rorth = zeros(size(R')) ;
% % for idim =1:3
% %     INDLOC =idim:3:size(R,1) ;
% %     Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ;
% % end
% 
% %
% % [Gb,dR,DOFr,DOFm] = BCs_BEAMS_PERIODIC(DOFA,DOFB,R,a_A,a_B,Rorth') ;
% 
% end