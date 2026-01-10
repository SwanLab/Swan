function [F,fextRVEr_glo,rRB_glo,fextDOMred_glo,FORCES_2_PRINT_glo]= ...
    AssemblyFdomains_vectorRVE_plate(DATAROM,MESH2D,DATAIN,FORCES,ndim)




% T1/



nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)
nmodes(1).REACTION=[] ;
nmodes(1).DISPLACEMENT=[] ;
for i = 1:length(DATAROM)
    
    nmodes(i).REACTION = size(DATAROM{i} .BasisRdef,2) ;
    nmodes(i).DISPLACEMENT = size(DATAROM{i} .BasisUdef,2) ;
end

% -----------------------
% LOOP OVER DOMAINS
% -----------------------
% Since all domains are assumed to be identical, the contribution to the
% external force vector due to gravity forces is computed just once


DATAIN = DefaultField(DATAIN,'INCLUDE_GRAVITY',1) ;
if DATAIN.INCLUDE_GRAVITY == 1
    FORCES = DefaultField(FORCES,'GRAVITY', [0 -9.81  0]) ;
    for istructure = 1:length(DATAROM)
        % Fgravity{istructure} = cell(nnodeE,1);
        % for inode = 1:nnodeE
        Fgravity{istructure} = DATAROM{istructure}.GAMMAforce*FORCES.GRAVITY(:) ;
        % end
        %   Fgravity{istructure} = cell2mat(Fgravity{istructure}) ;
        fextRVEr_gravity{istructure}  = DATAROM{istructure}.Lgravity*FORCES.GRAVITY(:)   ;
        rRB_gravity{istructure} = DATAROM{istructure}.BETA_gravity*FORCES.GRAVITY(:)   ;
        fextDOMred_gravity{istructure} =  DATAROM{istructure}.Dgravity*FORCES.GRAVITY(:)   ;
    end
else
    nstructures = length(DATAROM) ;
    Fgravity = cell(1,nstructures);
    Fgravity(:) = {0} ;
    fextRVEr_gravity  = Fgravity  ;
    rRB_gravity = Fgravity ;
    fextDOMred_gravity = Fgravity;
end



FORCES.TRACTIONS = DefaultField(FORCES.TRACTIONS,'ISLOCAL',zeros(size(FORCES.TRACTIONS.LOADS))) ;

%DATAIN = DefaultField(DATAIN,'PRINT_DISTRIBUTED_FORCES',0)  ;

FORCES_2_PRINT = cell(size(FORCES.TRACTIONS.LOADS)) ;  % Traction forces on each slice and each surface (per unit area)

nsurfaces = size(FORCES.TRACTIONS.LOADS,1) ;  % Each surface is assumed to have a uniform load
nsurfaces2 = length(MESH2D.ELEMENTS_SURFACES) ;
if nsurfaces2 == 0
    MESH2D.ELEMENTS_SURFACES{1} = [1:size(MESH2D.CN,1)]';
    nsurfaces2 = 1;
end
nsurfaces =  min(nsurfaces,nsurfaces2) ;



STRUCTURES = unique(MESH2D.MaterialType ) ;
nstructures = length(STRUCTURES) ;
F  =zeros(ndim*nnode,1) ;
fextRVEr = cell(nsurfaces,nstructures) ; % Variable R^*
fextDOMred =  cell(nsurfaces,nstructures) ; % Variable F*
rRB = cell(nsurfaces,nstructures) ; % Amplitude resultant reactions
Fbe = cell(nsurfaces,nstructures) ;

% Loop over lines
for e = 1:nsurfaces
    %    CNlocNOD = MESH2D.CN(e,:) ;
    %   CNloc = Nod2DOF(CNlocNOD,ndim) ;
    %
    % Computation of traction forces
    for ielemtype=1:length(STRUCTURES)
        elemtype = STRUCTURES(ielemtype) ;
        
        [FtractionE fextRVEr_traction,rRB_traction,fextDOMred_traction,FORCES_2_PRINT_local]= ...
            TractionForcesDomainRVE_plate(DATAROM{elemtype}.UPSILONforceTRAC,FORCES.TRACTIONS.LOADS(e,:),...
            FORCES.TRACTIONS.COOR(e,:),nnodeE,ndim,DATAROM{elemtype}.Q_tractions,...
            DATAROM{elemtype}.BETA_tractions,DATAROM{elemtype}.D_tractions,nmodes(elemtype),...
            FORCES.TRACTIONS.ISLOCAL(e,:)) ;
        nfacesLOC = min(size(FORCES_2_PRINT,2),length(FORCES_2_PRINT_local)) ;
        
        FORCES_2_PRINT(e,1:nfacesLOC) = FORCES_2_PRINT_local(1:nfacesLOC) ;
        fextRVEr{e,ielemtype} = fextRVEr_traction ; %+ fextRVEr_gravity{elemtype} ;
        fextDOMred{e,ielemtype} = fextDOMred_traction ; %+ fextDOMred_gravity{elemtype} ;
        Fbe{e,ielemtype} =   FtractionE ; % + Fgravity{elemtype};
        rRB{e,ielemtype} =   rRB_traction ; %rRB_gravity{elemtype} ;
    end
    
end















% Construction of elemental arrays
% --------------------------------
%F,fextRVEr,rRB,fextDOMred,FORCES_2_PRINT
fextDOMred_glo = cell(nelem,1) ;
ssCC = size(fextDOMred{1,1}) ;
fextDOMred_glo(:) = {zeros(ssCC)} ;

fextRVEr_glo = cell(nelem,1) ;
ssDD = size(fextRVEr{1,1}) ;
fextRVEr_glo(:) = {zeros(ssDD)} ;
Fbe_glo = cell(nelem,1) ;
ssAA = size(Fbe{1,1}) ;
Fbe_glo(:) = {zeros(ssAA)} ;
rRB_glo = cell(nelem,1) ;
ssBB= size(rRB{1,1}) ;
rRB_glo(:) = {zeros((ssBB))} ;

FORCES_2_PRINT_glo = cell(nelem,size(FORCES.TRACTIONS.LOADS,2)) ;
for isurface = 1:nsurfaces
    ELEM_SURFACE = MESH2D.ELEMENTS_SURFACES{isurface} ;
    for istructure =1:nstructures
        ELEMS = find(MESH2D.MaterialType(ELEM_SURFACE)==istructure); % Elements pertaining to this line and also to entity "istructure"
        ELEMS_LOC = ELEM_SURFACE(ELEMS) ;
        %  CUMULATIVE = 1 ;
        FA = {Fbe{isurface,istructure}}  ;
        FB = {rRB{isurface,istructure}}  ;
        FC = {fextDOMred{isurface,istructure}} ;
        FD =  {fextRVEr{isurface,istructure}}  ;
        
        if ~isempty(ELEMS_LOC)
            %             if CUMULATIVE == 0
            %                 Fbe_glo(ELEMS_LOC) =  {Fbe{isurface,istructure}}  ;
            %                 rRB_glo(ELEMS_LOC) = {rRB{isurface,istructure}}  ;
            %                 fextDOMred_glo(ELEMS_LOC) = {fextDOMred{isurface,istructure}}  ;
            %                 fextRVEr_glo(ELEMS_LOC) = {fextRVEr{isurface,istructure}}  ;
            %             else
            Fbe_glo = AddCellsValues(Fbe_glo,ELEMS_LOC,FA,ssAA) ;
            rRB_glo = AddCellsValues(rRB_glo,ELEMS_LOC,FB,ssBB) ;
            fextDOMred_glo = AddCellsValues(fextDOMred_glo,ELEMS_LOC,FC,ssCC) ;
            fextRVEr_glo = AddCellsValues(fextRVEr_glo,ELEMS_LOC,FD,ssDD) ;
            % end
        end
        
        
        
        for iface =1:size(FORCES.TRACTIONS.LOADS,2)
            FORCES_2_PRINT_glo(ELEMS_LOC,iface) = {FORCES_2_PRINT{isurface,istructure}}  ;
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADDING GRAVITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DATAIN.INCLUDE_GRAVITY == 1
    for istructure =1:nstructures
        %  elemtype = STRUCTURES(ielemtype) ;
        ELEMS = find(MESH2D.MaterialType ==istructure); % Elements pertaining to this line and also to entity "istructure"
        
        %  CUMULATIVE = 1 ;
        FA = { Fgravity{istructure}}  ;
        FB = {rRB_gravity{istructure}}  ;
        FC = {fextDOMred_gravity{istructure}} ;
        FD =  {fextRVEr_gravity{istructure}}  ;
        
        %            fextRVEr{e,ielemtype} = fextRVEr_traction ; %+ fextRVEr_gravity{elemtype} ;
        %         fextDOMred{e,ielemtype} = fextDOMred_traction ; %+ fextDOMred_gravity{elemtype} ;
        %         Fbe{e,ielemtype} =   FtractionE ; % + Fgravity{elemtype};
        %         rRB{e,ielemtype} =   rRB_traction ; %rRB_gravity{elemtype} ;
        
        %if ~isempty(ELEMS)
        
        Fbe_glo = AddCellsValues(Fbe_glo,ELEMS,FA,ssAA) ;
        rRB_glo = AddCellsValues(rRB_glo,ELEMS,FB,ssBB) ;
        fextDOMred_glo = AddCellsValues(fextDOMred_glo,ELEMS,FC,ssCC) ;
        fextRVEr_glo = AddCellsValues(fextRVEr_glo,ELEMS,FD,ssDD) ;
        
        %end
        
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% ASSEMBLY
% -----------------
%
%    for anod=1:nnodeE
%         a = Nod2DOF(anod,ndim) ;
%         Anod = MESH2D.CN(e,anod) ; A = Nod2DOF(Anod,ndim) ;
%         F(A) = F(A) + Fbe(a) ;
%     end

% EMPTYINDEX_1 = cellfun(@isempty,Fbe_glo) ;
% EMPTYINDEX = find(EMPTYINDEX_1==1) ;
% NEM = find(EMPTYINDEX_1==0) ;
%
% Fbe_glo(EMPTYINDEX) = {zeros(size(Fbe_glo{NEM(1)}))};
% rRB_glo(EMPTYINDEX) = {zeros(size(rRB_glo{NEM(1)}))};
% fextDOMred_glo(EMPTYINDEX) = {zeros(size(fextDOMred_glo{NEM(1)}))};
% fextRVEr_glo(EMPTYINDEX) = {zeros(size(fextRVEr_glo{NEM(1)}))};

Fbe = cell2mat(Fbe_glo) ;

for anod=1:nnodeE
    a = Nod2DOFelem(anod,ndim,nnodeE,nelem) ;
    Anod = MESH2D.CN(:,anod) ; A = Nod2DOF(Anod,ndim) ;
    F(A) = F(A) + Fbe(a) ;
end

end
% 
% function Fbe_glo = AddCellsValues(Fbe_glo,ELEMS_LOC,FA,ssAA)
% 
% FA = cell2mat(Fbe_glo(ELEMS_LOC)) + repmat(cell2mat(FA),length(ELEMS_LOC),1) ;
% FA = mat2cell(FA,ssAA(1)*ones(length(ELEMS_LOC),1),1) ;
% Fbe_glo(ELEMS_LOC) = FA ;
% 
% end
% 
% end
