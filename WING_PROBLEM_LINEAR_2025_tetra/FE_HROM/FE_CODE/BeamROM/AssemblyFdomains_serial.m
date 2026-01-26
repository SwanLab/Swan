function [F,fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    AssemblyFdomains_serial(DATAROM,MESH1D,DATAIN,FORCES_INPUT,ndim)

if nargin == 0
    load('tmp2.mat')
    DATAIN.INCLUDE_GRAVITY  = 1; 
end



% T1/
DATAIN = DefaultField(DATAIN,'TRACTION_FORCES_DEFINED_BY_LINES',0)
FORCES = FORCES_INPUT;

if DATAIN.TRACTION_FORCES_DEFINED_BY_LINES == 1
    % External forces defined line-wise, rather than element-wise
    nelem = size(MESH1D.CN,1) ;
    nfacesMAX = size(FORCES_INPUT.TRACTIONS.LOADS,2) ;
    FORCES.TRACTIONS.LOADS = cell(nelem,nfacesMAX) ; % Force per unit area for each slice and face
    FORCES.TRACTIONS.COOR = cell(nelem,nfacesMAX) ; % Coordinates of the points at which the forces
    FORCES.TRACTIONS.ISLOCAL = zeros(nelem,nfacesMAX) ; % Coordinates of the points at which the forces
    
    
    nlines = size(FORCES_INPUT.TRACTIONS.LOADS,1) ;
    for iline = 1:nlines
        % Nodes pertaining to this line
        NODES_LINE = MESH1D.NODES_LINES{iline} ;
        % Elements corresponding to these nodes
        [~,ELEM_LINE ]= ElemBnd(MESH1D.CN,NODES_LINE) ;
        for iface = 1:nfacesMAX
            if ~isempty(FORCES_INPUT.TRACTIONS.LOADS{iline,iface} )
            FORCES.TRACTIONS.LOADS(ELEM_LINE,iface) =  {FORCES_INPUT.TRACTIONS.LOADS{iline,iface}} ;
            end
            if ~isempty(FORCES_INPUT.TRACTIONS.LOADS{iline,iface} )
            FORCES.TRACTIONS.COOR(ELEM_LINE,iface) = { FORCES_INPUT.TRACTIONS.COOR{iline,iface}} ;
            end
                 
                FORCES.TRACTIONS.ISLOCAL(ELEM_LINE,iface) =  FORCES_INPUT.TRACTIONS.ISLOCAL(iline,iface) ;
             
        end
    end
    
end


nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)
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
        Fgravity{istructure} = cell(nnodeE,1);
        for inode = 1:nnodeE
            Fgravity{istructure}{inode} = DATAROM{istructure}.GAMMAforce{inode}*FORCES.GRAVITY(:) ;
        end
        Fgravity{istructure} = cell2mat(Fgravity{istructure}) ;
        fextBEAMr_gravity{istructure}  = DATAROM{istructure}.Lgravity*FORCES.GRAVITY(:)   ;
        rRB_gravity{istructure} = DATAROM{istructure}.BETA_gravity*FORCES.GRAVITY(:)   ;
        fextDOMred_gravity{istructure} =  DATAROM{istructure}.Dgravity*FORCES.GRAVITY(:)   ;
    end
else
    nstructures = length(DATAROM) ;
    Fgravity = cell(1,nstructures); 
    Fgravity(:) = {0} ; 
    fextBEAMr_gravity  = Fgravity  ;
    rRB_gravity = Fgravity ;
    fextDOMred_gravity = Fgravity;
end

F  =zeros(ndim*nnode,1) ;
fextBEAMr = cell(nelem,1) ; % Variable R^*
fextDOMred =  cell(nelem,1) ; % Variable F*
rRB = cell(nelem,1) ; % Amplitude resultant reactions

FORCES.TRACTIONS = DefaultField(FORCES.TRACTIONS,'ISLOCAL',zeros(size(FORCES.TRACTIONS.LOADS))) ;

%DATAIN = DefaultField(DATAIN,'PRINT_DISTRIBUTED_FORCES',0)  ;

FORCES_2_PRINT = cell(size(FORCES.TRACTIONS.LOADS)) ;  % Traction forces on each slice and each surface (per unit area)


for e = 1:nelem
    %    CNlocNOD = MESH1D.CN(e,:) ;
    %   CNloc = Nod2DOF(CNlocNOD,ndim) ;
    %
    % Computation of traction forces
    elemtype = MESH1D.MaterialType(e) ;
    
    [FtractionE, fextBEAMr_traction,rRB_traction,fextDOMred_traction,FORCES_2_PRINT_local]= ...
        TractionForcesDomain(DATAROM{elemtype}.UPSILONforceTRAC,FORCES.TRACTIONS.LOADS(e,:),...
        FORCES.TRACTIONS.COOR(e,:),nnodeE,ndim,DATAROM{elemtype}.Q_tractions,...
        DATAROM{elemtype}.BETA_tractions,DATAROM{elemtype}.D_tractions,nmodes(elemtype),...
        FORCES.TRACTIONS.ISLOCAL(e,:)) ;
    nfacesLOC = min(size(FORCES_2_PRINT,2),length(FORCES_2_PRINT_local)) ;
    FORCES_2_PRINT(e,1:nfacesLOC) = FORCES_2_PRINT_local(1:nfacesLOC) ;
    fextBEAMr{e} = fextBEAMr_traction + fextBEAMr_gravity{elemtype} ;
    fextDOMred{e} = fextDOMred_traction + fextDOMred_gravity{elemtype} ;
    Fbe =  Fgravity{elemtype} + FtractionE ;
    rRB{e} = rRB_gravity{elemtype} + rRB_traction ;
    for anod=1:nnodeE
        a = Nod2DOF(anod,ndim) ;
        Anod = MESH1D.CN(e,anod) ; A = Nod2DOF(Anod,ndim) ;
        F(A) = F(A) + Fbe(a) ;
    end
end