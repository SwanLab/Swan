function [qDEF,qRB]= AmplitudeDisplacements_RVE(DATAROM,MESH2D,rDEF,fextDOMred,...
    DATA_REFMESH,a,DATAIN,ndim,DOFsKEEP) ;


if nargin ==0
    load('tmp3.mat')
end

nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)

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
nMODESrb = size(DATAROM{1}.BasisIntRB{1},2) ;

qRB = zeros(nMODESrb,size(qDEF,2)) ;
for ientity = 1:nentities
    ELEMS = find(MESH2D.MaterialType == ientity) ;
    
    
    f = DATAROM{ientity}.f ;
    fI= DATAROM{ientity}.fI ;
    R_f = DATA_REFMESH{ientity}.BasisUrb(f,:) ; % Rigid body modes faces 1 and 2
    M = DATA_REFMESH{ientity}.M ;
    
    
    
    if ~isstruct(DATAROM{ientity}.BasisInt)
        V = DATAROM{ientity}.BasisInt ;
        T = cell(size(V)) ;
        for iface = 1:nnodeE
            
            
            RotationLocal = [] ; 
            if isfield(DATA_REFMESH{ientity},'RotationMatrixFace')
                RotationLocal =  DATA_REFMESH{ientity}.RotationMatrixFace{iface}' ; 
            end
            
            BasisRrb_local = M(fI{iface},fI{iface})*DATA_REFMESH{ientity}.BasisUrb(fI{iface},:) ;
         %   RotationLocal = ;
            T{iface}  = RotateMatricesProduct_interface(RotationLocal,BasisRrb_local,V{iface}) ;
            
            %  T{iface} = V{iface}'*M(fI{iface},fI{iface})*DATA_REFMESH{ientity}.BasisUrb(fI{iface},:) ;
            %   T_2 = V2'*M(f2,f2)*DATA_REFMESH{ientity}.BasisUrb(f2,:) ;
        end
        
        
        RRinv = inv(R_f'*M(f,f)*R_f) ;
        C_q = -RRinv*(R_f'*M(f,f)*DATAROM{ientity}.BasisUdef(f,:)) ;
        C_a = cell(size(V)) ;
        for iface = 1:nnodeE
            C_a{iface} = RRinv*T{iface}' ;  % Corrected 16-May-2018-
            % C_a{2} = RRinv*T_2' ;
        end
        
    else
        
        % New method, kinematically constrained
        V = DATAROM{ientity}.BasisInt.BasisINTF ;
        [DOFsNOD,~]=  cellfun(@size,DATAROM{ientity}.BasisInt.BasisINTFall_cell) ;
        T  = V'*M(f,f)*DATA_REFMESH{ientity}.BasisUrb(f,:) ;
        RRinv = inv(R_f'*M(f,f)*R_f) ;
        C_q = -RRinv*(R_f'*M(f,f)*DATAROM{ientity}.BasisUdef(f,:)) ;
        C_a = RRinv*T' ;  %
        
    end
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qRB(:,ELEMS) = C_q*qDEF(:,ELEMS);
    
    CNlocNOD = MESH2D.CN(ELEMS,:) ;
    
    
    if ~isstruct(DATAROM{ientity}.BasisInt)
        
        % old method
        
        if any(ndim - ndim(1))
            aALL = zeros(nnode*max(ndim),1) ;
            aALL(DOFsKEEP) = a;
            for inode = 1:nnodeE
                NODE = CNlocNOD(:,inode) ;
                DOFS = Nod2DOF(NODE,max(ndim)) ;
                aDOFS = aALL(DOFS) ;
                aDOFS = reshape(aDOFS,max(ndim),[]) ;
                %  rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{inode,ientity}*aDOFS  ;
                nrows = size( C_a{inode},1) ;
                C_aLOC = zeros(nrows,max(ndim)) ;
                C_aLOC(1:nrows,1:ndim(inode)) =  C_a{inode} ;
                qRB(:,ELEMS) = qRB(:,ELEMS) +C_aLOC*aDOFS  ;
            end
        else
            ndim = ndim(1) ;
            for inode = 1:nnodeE
                NODE = CNlocNOD(:,inode) ;
                DOFS = Nod2DOF(NODE,ndim) ;
                aDOFS = a(DOFS) ;
                aDOFS = reshape(aDOFS,ndim,[]) ;
                %  rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{inode,ientity}*aDOFS  ;
                qRB(:,ELEMS) = qRB(:,ELEMS) + C_a{inode}*aDOFS  ;
            end
            
        end
        
    else
         % New method, Apr-2019. Kinematically constrained
        
        if any(ndim - ndim(1))
            aALL = zeros(nnode*max(ndim),1) ;
            aALL(DOFsKEEP) = a;
            COLUMNS_ini = 1; 
            for inode = 1:nnodeE
                NODE = CNlocNOD(:,inode) ;
                DOFS = Nod2DOF(NODE,max(ndim)) ;
                aDOFS = aALL(DOFS) ;
                aDOFS = reshape(aDOFS,max(ndim),[]) ;
                %  rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{inode,ientity}*aDOFS  ;
                nrows = size( C_a,1) ;
                C_aLOC = zeros(nrows,max(ndim)) ;
                COLUMNS_fin = COLUMNS_ini +  ndim(inode)-1 ;
                COLUMNS = COLUMNS_ini:COLUMNS_fin ;
                COLUMNS_ini = COLUMNS_fin + 1 ; 
                C_aLOC(1:nrows,1:ndim(inode)) =  C_a(:,COLUMNS)  ;
                
                 
                
                qRB(:,ELEMS) = qRB(:,ELEMS) +C_aLOC*aDOFS  ;
            end
        else
            ndim = ndim(1) ;
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
            
        end
        
    end
    
    
    
    
end


end
