function [K,DOFsKEEP] = AssemblyKrve_notequal(DATAROM,MESH2D,DATAIN,ndim)
% JAHO, 19-Nov-2018
% Reduced-order Stiffness matrix
% Interfaces may have different number of modes

if nargin == 0
    load('tmp2.mat')
end

%% We have that ndim(1) = ndim(3), and ndim(2) = 4
%  The idea is to make ndimMAX = max(ndim), and proceed as in the normal
%  case, but filling Kskel with zeros (so that it is 4 ndimMAX x 4 ndimMAX)  
% ------------------------------------------------------------------
ndimMAX = max(ndim) ;
[DOFSeliminate ]= DOFSeliminateNEQDIM(ndim,MESH2D) ; 

 
% ---------------------------------------------
for itype = 1:length(DATAROM)
    Ke = DATAROM{itype}.Kskel ;
    Ke = mat2cell(Ke,ndim,ndim); 
    for iface = 1:size(Ke,1)
        for jface =1:size(Ke,2) 
            [nrows,ncols] = size(Ke{iface,jface}) ;  
            KeLOCAL = zeros(ndimMAX,ndimMAX) ;
            KeLOCAL(1:nrows,1:ncols) = Ke{iface,jface} ; 
            Ke{iface,jface} = KeLOCAL ; 
        end
    end
     DATAROM{itype}.Kskel = cell2mat(Ke);    
end

 

nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
%Ke = DATAROM.Kskel ; % "Elemental "Stiffness matrix  (the same for all domains)
nnodeE = length(ndim); % Number of nodes per element (number of interfaces per element)

DATAIN = DefaultField(DATAIN,'VECTORIZED_CODE',1) ;  
ndim = ndimMAX ; 
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


% Now we remove the indices DOFSeliminate 
% -----------------------------------------
DOFSeliminate = sort(DOFSeliminate) ; 
DOFsKEEP = setdiff(1:size(K,1),DOFSeliminate) ; 
K = K(DOFsKEEP,DOFsKEEP) ; 



end