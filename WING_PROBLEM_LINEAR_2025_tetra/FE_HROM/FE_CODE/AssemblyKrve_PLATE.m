function [K,ndim] = AssemblyKrve_PLATE(DATAROM,MESH2D,DATAIN,DATA_REFMESH)
% JAHO, 6-September-2018 
% Reduced-order Stiffness matrix. Plate structures
% All interfaces are assumed to have the same number of modes 

if nargin == 0
    load('tmp2.mat')
end

 
nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
nnodeE = size(MESH2D.CN,2); % Number of nodes per element (number of interfaces per element)
ndim = size(DATA_REFMESH{1}.BasisINTfc,2)/nnodeE  ;


DATAIN = DefaultField(DATAIN,'VECTORIZED_CODE',1) ;  

if DATAIN.VECTORIZED_CODE == 0
    
    % Non-vectorized version
    K = sparse(nnode*ndim,nnode*ndim) ;

    for e = 1:nelem
        elemtype = MESH2D.MaterialType(e) ; 
        Ke = DATAROM{elemtype}.Kskel ; 
        for anod=1:nnodeE
            a = Nod2DOF(anod,ndim) ;
            for bnod= 1:nnodeE
                b = Nod2DOF(bnod,ndim) ;
                Anod = MESH2D.CN(e,anod) ;  A = Nod2DOF(Anod,ndim) ;
                Bnod = MESH2D.CN(e,bnod) ;  B = Nod2DOF(Bnod,ndim) ;
                %%%%%
                K(A,B) = K(A,B) + Ke(a,b) ;
            end
        end
    end
else
    % Vectorized version 
    % ------------------
    Kelem = cell(nelem,1) ; 
    % Stiffness matrix containing all local stiffness
    for itype = 1:length(DATAROM)
         Ke = DATAROM{itype}.Kskel ; 
         ELEMS = find(MESH2D.MaterialType == itype) ; 
         Kelem(ELEMS) = {Ke} ; 
          
    end
    Kelem = cell2mat(Kelem) ; 
    
    
    
    
    %     
    m = nnode*ndim ; % Number of rows
    n = m ;          % Number of columns
    nzmaxLOC = size(Kelem,1)*size(Kelem,2) ;   % Maximum number of zeros (number of entries of Belems)
    K = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Bst
    
    for anod=1:nnodeE % Loop over element nodes (rows)
        a = Nod2DOFelem(anod,ndim,nnodeE,nelem) ;  % ROWS number (in Kelem) for node   "anod" (for all elements)
        for bnod= 1:nnodeE  % Loop over element nodes (columns)
            b = Nod2DOF(bnod,ndim) ;
            Anod = MESH2D.CN(:,anod) ;  A = Nod2DOF(Anod,ndim) ;  % DOFs in the global K matrix
            Bnod = MESH2D.CN(:,bnod) ;  B = Nod2DOF(Bnod,ndim) ;
            %%%%%
            %  K(A,B) = K(A,B) + Kelem(a,b) ;
            %%%%%
            s = Kelem(a,b) ;
            s=s(:) ;  nzmax = length(s);
            % Indices "i" and "j"
            i = repmat(A,ndim,1);
            j = ObtainIndexJassemb(B,ndim) ;
            %%%%
            
            K = K + sparse(i,j,s,m,n,length(s)) ;
        end
    end
    
    
end


end