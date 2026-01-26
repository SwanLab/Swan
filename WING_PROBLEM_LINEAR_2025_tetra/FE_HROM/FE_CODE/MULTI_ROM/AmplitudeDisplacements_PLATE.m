function [qDEF,qRB]= AmplitudeDisplacements_PLATE(DATAROM,MESH2D,rDEF,fextDOMred,DATA_REFMESH,a,DATAIN,ndim) ;


if nargin ==0
    load('tmp0.mat')
end

nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
nnodeE = size(MESH2D.CN,2); % Number of nodes per element (number of interfaces per element)

nentities = length(DATAROM) ;

nmodesU = size(fextDOMred{1},1) ; % Numbmer of displacement modes
nmodesR = size(rDEF,1); % Number of reaction modes

F = cell2mat(fextDOMred) ;  % Reduced external froces
F = reshape(F,nmodesU,[]) ;

% rDEF = cell2mat(rDEF) ;
% rDEF = reshape(rDEF,nmodesR,[]) ;  % Amplitude self-equilibrated reactions
qDEF = zeros(size(F)) ;

for ientity = 1:nentities
    ELEMS = find(MESH2D.MaterialType == ientity) ;
    KHinv = DATAROM{ientity}.KdomRED\DATAROM{ientity}.Hqr ;
    qDEF(:,ELEMS) = KHinv*rDEF(:,ELEMS)  +DATAROM{ientity}.KdomRED\F(:,ELEMS) ; % Amplitude deformational displacements
    
end



%%% RIGID BODY AMPLITUDES
nMODESrb = 6 ;

qRB = zeros(nMODESrb,size(qDEF,2)) ;
for ientity = 1:nentities
    ELEMS = find(MESH2D.MaterialType == ientity) ;
    
    
    %f1 = DATAROM{ientity}.f1 ;  % DOFs face 1
    %f2 = DATAROM{ientity}.f2 ;  % DOFs face 2
    f = DATAROM{ientity}.f ;
   % fI= DATAROM{ientity}.fI ;
    R_f = DATA_REFMESH{ientity}.BasisUrb(f,:) ; % Rigid body modes faces 1 and 2
    M = DATA_REFMESH{ientity}.M ;
    
    %     if iscell(DATAROM{ientity}.BasisInt)
    %         V1 = DATAROM{ientity}.BasisInt{1} ;
    %         V2 = DATAROM{ientity}.BasisInt{2} ;
    %     else
    %         V1 = DATAROM{ientity}.BasisInt ;
    %         V2 = V1 ;
    %     end
   % V = DATAROM{ientity}.BasisInt ;
    V = DATA_REFMESH{ientity}.BasisINTfc ;  
 %   T = cell(size(V)) ;
  %  for iface = 1:nnodeE
        T  = V'*M(f,f)*DATA_REFMESH{ientity}.BasisUrb(f,:) ;
        %   T_2 = V2'*M(f2,f2)*DATA_REFMESH{ientity}.BasisUrb(f2,:) ;
   % end
    
    RRinv = inv(R_f'*M(f,f)*R_f) ;
    C_q = -RRinv*(R_f'*M(f,f)*DATAROM{ientity}.BasisUdef(f,:)) ;
   % C_a = cell(size(V)) ;
   % for iface = 1:nnodeE
        C_a = RRinv*T' ;  % 
        % C_a{2} = RRinv*T_2' ;
   % end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qRB(:,ELEMS) = C_q*qDEF(:,ELEMS);
    
    CNlocNOD = MESH2D.CN(ELEMS,:) ;
    for inode = 1:nnodeE
        NODE = CNlocNOD(:,inode) ;
        DOFS = Nod2DOF(NODE,ndim) ;
        aDOFS = a(DOFS) ;
        aDOFS = reshape(aDOFS,ndim,[]) ;
        %  rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{inode,ientity}*aDOFS  ;
         COLUMNS_ini =  (inode-1)*ndim+1 ; 
           COLUMNS_fin = inode*ndim ;  
         COLUMNS = COLUMNS_ini:COLUMNS_fin ; 
        qRB(:,ELEMS) = qRB(:,ELEMS) + C_a(:,COLUMNS)*aDOFS  ;
    end
    
    %     for ielem = 1:nelem
    %         CNlocNOD = MESH2D.CN(ielem,:) ;
    %         % CNloc = Nod2DOF(CNlocNOD,ndim) ;
    %         qRB(:,ielem) = C_q*qDEF(:,ielem);
    %         for inode = 1:nnodeE
    %             NODE = CNlocNOD(inode) ;
    %             DOFS = Nod2DOF(NODE,ndim) ;
    %             qRB(:,ielem) = qRB(:,ielem) + C_a{inode}*a(DOFS)  ;
    %         end
    %     end
    
    
    
    %
    %     nmodes = size(fextBEAMr{1},1) ;   % number of modes reactions
    % fextBEAMr = cell2mat(fextBEAMr)  ;
    % rDEF = -fextBEAMr ;
    % rDEF  =reshape(rDEF,nmodes,[]) ;
    % for  ientity = 1:length(DATAROM)
    %     ELEMS = find(MESH2D.MaterialType == ientity) ;
    %     CNlocNOD = MESH2D.CN(ELEMS,:) ;
    %     for inode = 1:nnodeE
    %         NODE = CNlocNOD(:,inode) ;
    %         DOFS = Nod2DOF(NODE,ndim) ;
    %         aDOFS = a(DOFS) ;
    %         aDOFS = reshape(aDOFS,ndim,[]) ;
    %         rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{inode,ientity}*aDOFS  ;
    %     end
    %
    % end
    
    
    
end



% qDEF = qDEFmatr;
% qRB  =cell2mat(qRB);



%
%
% DATAIN = DefaultField(DATAIN,'RECOMPUTE_RIGID_BODY_DISPLACEMENTS',0) ;
%
% if DATAIN.RECOMPUTE_RIGID_BODY_DISPLACEMENTS == 1
%     qRB =    RecomputeRigidBodyAmplitudes(qDEF,DATA_REFMESH,DATAROM,a) ;
% end
end
