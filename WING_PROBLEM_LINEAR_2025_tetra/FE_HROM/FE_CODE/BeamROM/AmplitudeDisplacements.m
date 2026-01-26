function [qDEF,qRB]= ...
    AmplitudeDisplacements(DATAROM,MESH1D,rDEF,fextDOMred,DATA_REFMESH,a,DATAIN,ndim) ;


if nargin ==0
    load('tmp1.mat')
end

nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)

nentities = length(DATAROM) ;

nmodesU = size(fextDOMred{1},1) ; % Numbmer of displacement modes
nmodesR = size(rDEF,1); % Number of reaction modes

F = cell2mat(fextDOMred) ;  % Reduced external froces
F = reshape(F,nmodesU,[]) ;

% rDEF = cell2mat(rDEF) ;
% rDEF = reshape(rDEF,nmodesR,[]) ;  % Amplitude self-equilibrated reactions
qDEF = zeros(size(F)) ;

for ientity = 1:nentities
    ELEMS = find(MESH1D.MaterialType == ientity) ;
    KHinv = DATAROM{ientity}.KdomRED\DATAROM{ientity}.Hqr ;
    qDEF(:,ELEMS) = KHinv*rDEF(:,ELEMS)  +DATAROM{ientity}.KdomRED\F(:,ELEMS) ; % Amplitude deformational displacements
    
end



%%% RIGID BODY AMPLITUDES
nMODESrb = size(DATAROM{1}.BasisIntRB{1},2);

qRB = zeros(nMODESrb,size(qDEF,2)) ;
for ientity = 1:nentities
    ELEMS = find(MESH1D.MaterialType == ientity) ;
    
    f1 = DATAROM{ientity}.f1 ;  % DOFs face 1
    f2 = DATAROM{ientity}.f2 ;  % DOFs face 2
    f = [f1;f2] ;
    R_f = DATA_REFMESH{ientity}.BasisUrb(f,:) ; % Rigid body modes faces 1 and 2
    M = DATA_REFMESH{ientity}.M ;
    
    if iscell(DATAROM{ientity}.BasisInt)
        V1 = DATAROM{ientity}.BasisInt{1} ;
        V2 = DATAROM{ientity}.BasisInt{2} ;
    else
        V1 = DATAROM{ientity}.BasisInt ;
        V2 = V1 ;
    end
    
    %
    
    %
    % iface= 1;
    T_1 = V1'*M(f1,f1)*DATA_REFMESH{ientity}.BasisUrb(f1,:) ;
    % iface= 2;
    %  VB2 = V'*BasisRdef(f2,:) ;
    BasisRrb_f2 = M(f2,f2)*DATA_REFMESH{ientity}.BasisUrb(f2,:) ;
    DATA_REFMESH{ientity} = DefaultField(DATA_REFMESH{ientity},'RotationMatrixFace',[] ) ;
    if ~isempty(DATA_REFMESH{ientity}.RotationMatrixFace)
        RotationLocal = DATA_REFMESH{ientity}.RotationMatrixFace{2}' ;
        T_2 = RotateMatricesProduct_interface(RotationLocal,BasisRrb_f2,V2) ;
    else
        T_2 = V2'*M(f2,f2)*DATA_REFMESH{ientity}.BasisUrb(f2,:) ;
        
    end
    
    
    
    RRinv = inv(R_f'*M(f,f)*R_f) ;
    C_q = -RRinv*(R_f'*M(f,f)*DATAROM{ientity}.BasisUdef(f,:)) ;
    C_a{1} = RRinv*T_1' ;  % Corrected 16-May-2018-
    C_a{2} = RRinv*T_2' ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qRB(:,ELEMS) = C_q*qDEF(:,ELEMS);
    
    CNlocNOD = MESH1D.CN(ELEMS,:) ;
    for inode = 1:nnodeE
        NODE = CNlocNOD(:,inode) ;
        DOFS = Nod2DOF(NODE,ndim) ;
        aDOFS = a(DOFS) ;
        aDOFS = reshape(aDOFS,ndim,[]) ;
        %  rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{inode,ientity}*aDOFS  ;
        qRB(:,ELEMS) = qRB(:,ELEMS) + C_a{inode}*aDOFS  ;
    end
    
end

%%%% 




end
