function [Tdef,Trb,R] = ComputeOperatorsGeneralizedForces(MESH1D,DATAROM,ientity,nmodes,DATA_REFMESH,ELEMS)

%  % Nodes corresponding to face "1"
    M = DATA_REFMESH{ientity}.M ;
    
    if iscell(DATAROM{ientity}.BasisInt)
        V1 = DATAROM{ientity}.BasisInt{1} ;        
        V2 = DATAROM{ientity}.BasisInt{2} ;
    else
        V1 = DATAROM{ientity}.BasisInt ;
        V2 = V1;        
    end
    V2 = V2(:,1:nmodes) ; 
    V1 = V1(:,1:nmodes) ;
    f1 = DATAROM{ientity}.f1 ;
    f2 = DATAROM{ientity}.f2 ;
    Tdef{1} = V1'*DATAROM{ientity}.BasisRdef(f1,:) ;
    Trb{1} = V1'*M(f1,f1)*DATA_REFMESH{ientity}.BasisUrb(f1,:) ;
    
    % Face 2. Rotation matrix 
    % ------------------------
    DATA_REFMESH{ientity}  =DefaultField(DATA_REFMESH{ientity},'RotationMatrixFace',[]) ; 
    if  ~isempty(DATA_REFMESH{ientity}.RotationMatrixFace{2})
        R = DATA_REFMESH{ientity}.RotationMatrixFace{2};
    else
        R = []  ; 
    end
   % ndimSP = size(R,1) ; 
    % BasisUrb and BasisRdef are given in the domain reference system
    BasisRdef_2 = DATAROM{ientity}.BasisRdef(f2,:) ; 
    BasisRrb_2 = M(f2,f2)*DATA_REFMESH{ientity}.BasisUrb(f2,:) ;
    % Transformation to the reference system attached to face 2 
    BasisRdef_2 = RotateMatrix(R',BasisRdef_2) ; 
    BasisRrb_2 = RotateMatrix(R',BasisRrb_2) ; 

%     
    Tdef{2} = V2'*BasisRdef_2 ;
    Trb{2} = V2'*BasisRrb_2 ;