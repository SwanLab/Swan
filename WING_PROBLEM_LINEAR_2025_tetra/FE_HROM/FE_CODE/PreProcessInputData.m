function [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,celasglo,...
    DOFr,dR,Tnod,CNb,fNOD,Fpnt,NameFileMesh,densglo,celasgloINV] = ...
    PreProcessInputData(NameFileMeshDATA,PROPMAT,DIRICHLET,NEUMANN,POINT_FORCE,...
    fBODY,dens0,typePROBLEM)



% ----------------------------
%%%%%%%%%%%%%%%%%
 % ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB
NameFileMesh = [NameFileMeshDATA,'.msh']; % Name of the file containing the mesh information (Generated with GID)
[COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
    ReadMeshFile(NameFileMesh)  ;

nnode = size(COOR,1) ;% Number of nodes 
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. MATERIAL PROPERTIES: output celasglo   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%typePROBLEM = 'pstress';  %'pstress'/'pstrain'/'3D';  Plane stress/ plane strain problem 
if ndim==2 
    nstrain = 3; 
else
    nstrain = 6 ; 
    typePROBLEM ='3D' ;
end
celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)    
    celas3D =PROPMAT(imat).ElasticityMatrix ; %
    INVcelas3D = inv(celas3D) ; 
    ELEMS = find(MaterialType == imat) ;
    
    switch typePROBLEM
        case 'pstrain'
            rowcol = [1 2 6] ;
            celas = celas3D(rowcol,rowcol) ;
        case 'pstress'
            rowcol = [1 2 6] ;
            celasINV3D = inv(celas3D) ;
            celasINV = celasINV3D(rowcol,rowcol) ;
            celas = inv(celasINV) ;
        case '3D'
            celas = celas3D ;
    end
    for eLOC=1:length(ELEMS)
        e = ELEMS(eLOC) ;
        celasglo(:,:,e) = celas ;
        celasgloINV(:,:,e) = INVcelas3D ; 
    end
end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of nodes at which displacement is prescribed (in any of the x-y and z directions)
rnod = cell(ndim,1) ; uPRES=cell(ndim,1)  ; 
for  icond = 1:length(DIRICHLET)
    rnodLOC=  ListOfNodesFACES(NameFileMesh,DIRICHLET(icond).NUMBER_SURFACE,ndim) ;
    PRESCRIBED_DISPLACEMENT = DIRICHLET(icond).PRESCRIBED_DISP; 
    for idim = 1:ndim
        if ~isempty(PRESCRIBED_DISPLACEMENT{idim})
            displ1 = PRESCRIBED_DISPLACEMENT{idim} ; 
            rnod{idim} = [rnod{idim}; rnodLOC] ; 
            uPRES{idim} =  [uPRES{idim}; displ1*ones(size(rnod{idim}))] ;  
        end
    end       
end
% Removed repeated condions 
% ---------------------------
for idim = 1:ndim 
   [rnod{idim}, AAA] = unique(rnod{idim}) ;
   uPRES{idim} = uPRES{idim}(AAA) ; 
end
 
% Degrees of freedom and prescribed displacements 
DOFr = [] ; dR = [] ; 
for idim = 1:ndim 
    DOFr = [DOFr ; (rnod{idim}-1)*ndim+idim]; 
    dR = [dR ; uPRES{idim}]; 
end
 
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: Tnod, CNb, Fnod  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS
% ------------------------
CNb =cell(ndim,1) ; Tnod=cell(ndim,1) ;
for  icond = 1:length(NEUMANN)
    rnodLOC=  ListOfNodesFACES(NameFileMesh,NEUMANN(icond).NUMBER_SURFACE,ndim) ;
    FORCE_PER_UNIT_SURFACE = NEUMANN(icond).FORCE_PER_UNIT_SURFACE;
    for idim = 1:ndim
        if FORCE_PER_UNIT_SURFACE(idim) ~= 0
            t = FORCE_PER_UNIT_SURFACE(idim) ;
            CNbLOC = ElemBnd(CONNECTb,rnodLOC) ;
            TnodLOC = t*ones(size(CNbLOC)) ;
            CNb{idim} = [CNb{idim} ; CNbLOC] ;
            Tnod{idim} = [Tnod{idim}  ; TnodLOC] ;
        end
    end
end
      

% POINT LOADS 
% -----------
Fpnt =zeros(ndim*nnode,1) ; % There are no point loads
for iforce = 1:length(POINT_FORCE) 
    FORCE = POINT_FORCE(iforce).VALUE(1:ndim) ; 
    NODELOC = POINT_FORCE(iforce).NODE ; 
    DOFSloc =    Nod2DOF(NODELOC,ndim) ; 
    Fpnt(DOFSloc) = Fpnt(DOFSloc) + FORCE(:) ;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Body forces
%% 
fNOD = fBODY*ones(nnode*ndim,1) ; 

% 6. Density 
densglo = dens0*ones(nelem,1);  % in kKg/m^3 


% 
