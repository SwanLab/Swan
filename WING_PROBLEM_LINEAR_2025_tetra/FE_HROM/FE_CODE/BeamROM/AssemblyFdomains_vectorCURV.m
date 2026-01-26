function [F,fextBEAMr_glo,rRB_glo,fextDOMred_glo,FORCES_2_PRINT_glo]= ...
    AssemblyFdomains_vectorCURV(DATAROM,MESH1D,DATAIN,FORCES,ndim,DATA_REFMESH_glo)

if nargin == 0
    load('tmp4.mat')
end
FORCES_2_PRINT_glo = [];
ndimSP = size(MESH1D.COOR,2) ;


nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)
nmodes(1).REACTION=[] ;
nmodes(1).DISPLACEMENT=[] ;
for i = 1:length(DATAROM)
    
    nmodes(i).REACTION = size(DATAROM{i} .BasisRdef,2) ;
    nmodes(i).DISPLACEMENT = size(DATAROM{i} .BasisUdef,2) ;
end
STRUCTURES = unique(MESH1D.MaterialType ) ;

% ROTATIONS
% ---------

ROTATIONS_DOMAIN = mat2cell(MESH1D.ROTATIONS,ndimSP,ndimSP*ones(1,nelem)) ;
ROTATIONS_DOMAIN_transpose =  cellfun(@transpose,ROTATIONS_DOMAIN,'UniformOutput',false) ;  % Transpose


% -----------------------
% LOOP OVER DOMAINS
% -----------------------
% GRAVITY FORCES. 23-MAY-2019 
% ---------------------------------------------------------------
 [fextDOMred_glo,rRB_glo,fextBEAMr_glo,Fbe] = ...
     GravityForcesCurved(DATAIN,FORCES,ndimSP,nelem,ROTATIONS_DOMAIN_transpose,...
    DATAROM,MESH1D,nnodeE) ; 


 %%% TRACTION FORCES 
 % -----------------------
FORCES.TRACTIONS = DefaultField(FORCES.TRACTIONS,'ISLOCAL',zeros(size(FORCES.TRACTIONS.LOADS))) ;
FORCES_2_PRINT = cell(size(FORCES.TRACTIONS.LOADS)) ;  % No longer available
nlines = size(FORCES.TRACTIONS.LOADS,1) ;  % Each line is assumed to have a uniform load
nlinesGID = length(MESH1D.NODES_LINES) ; 
nlines = min(nlines,nlinesGID); 
nstructures = length(STRUCTURES) ;
F  =zeros(ndim*nnode,1) ;

% Loop over lines
istructure = 1; 

for e = 1:nlines    %
    NODES_LINE = MESH1D.NODES_LINES{e} ;
    [~,ELEM_LINE ]= ElemBnd(MESH1D.CN,NODES_LINE) ;
    
    for ielemtype=1:length(DATAROM)
        ELEMS = find(MESH1D.MaterialType(ELEM_LINE)==istructure); % Elements pertaining to this line and also to entity "istructure"
        ELEMS_LOC = ELEM_LINE(ELEMS) ;
        
        [FtractionE fextBEAMr_traction,rRB_traction,fextDOMred_traction,FORCES_2_PRINT_local]= ...
            TractionForcesDomain_curve(DATAROM{ielemtype}.UPSILONforceTRAC,FORCES.TRACTIONS.LOADS(e,:),...
            FORCES.TRACTIONS.COOR(e,:),nnodeE,ndim,DATAROM{ielemtype}.Q_tractions,...
            DATAROM{ielemtype}.BETA_tractions,DATAROM{ielemtype}.D_tractions,nmodes(ielemtype),...
            FORCES.TRACTIONS.ISLOCAL(e,:),ndimSP,ROTATIONS_DOMAIN(ELEMS_LOC),ROTATIONS_DOMAIN_transpose(ELEMS_LOC)) ;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Updating fextDOMred_glo
        fextDOMred_glo = UpdateCellForces(fextDOMred_glo,ELEMS_LOC,fextDOMred_traction) ; 
         % Updating Fbe
        Fbe = UpdateCellForces(Fbe,ELEMS_LOC,FtractionE) ; 
        % Updating rRB
        rRB_glo = UpdateCellForces(rRB_glo,ELEMS_LOC,rRB_traction) ; 
        % Updating fextBEAMr
        fextBEAMr_glo = UpdateCellForces(fextBEAMr_glo,ELEMS_LOC,fextBEAMr_traction) ; 
        %%%%%%%%%%%%%%%%%%%%%%
         
        
        
    end
    
end




% ASSEMBLY
fextBEAMr_glo = fextBEAMr_glo(:);
rRB_glo = rRB_glo(:);
fextDOMred_glo = fextDOMred_glo(:);
Fbe = cell2mat(Fbe(:));

for anod=1:nnodeE
    a = Nod2DOFelem(anod,ndim,nnodeE,nelem) ;
    Anod = MESH1D.CN(:,anod) ; A = Nod2DOF(Anod,ndim) ;
    F(A) = F(A) + Fbe(a) ;
end









%
% % --------------------------------
% % Construction of elemental arrays
% % --------------------------------
% %F,fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT
% fextDOMred_glo = cell(nelem,1) ;
% % ssCC = size(fextDOMred{1,1}) ;
% % fextDOMred_glo(:) = {zeros(ssCC)} ;
%
% fextBEAMr_glo = cell(nelem,1) ;
% ssDD = size(fextBEAMr{1,1}) ;
% fextBEAMr_glo(:) = {zeros(ssDD)} ;
% Fbe_glo = cell(nelem,1) ;
% ssAA = size(Fbe{1,1}) ;
% Fbe_glo(:) = {zeros(ssAA)} ;
% rRB_glo = cell(nelem,1) ;
% ssBB= size(rRB{1,1}) ;
% rRB_glo(:) = {zeros((ssBB))} ;
%
% FORCES_2_PRINT_glo = cell(nelem,size(FORCES.TRACTIONS.LOADS,2)) ;
% for iline = 1:nlines
%     NODES_LINE = MESH1D.NODES_LINES{iline} ;
%     %         % Elements corresponding to these nodes
%     [~,ELEM_LINE ]= ElemBnd(MESH1D.CN,NODES_LINE) ;
%     for istructure =1:nstructures
%         ELEMS = find(MESH1D.MaterialType(ELEM_LINE)==istructure); % Elements pertaining to this line and also to entity "istructure"
%         ELEMS_LOC = ELEM_LINE(ELEMS) ;
%         %  CUMULATIVE = 1 ;
%         FA = {Fbe{iline,istructure}}  ;
%         FB = {rRB{iline,istructure}}  ;
%         FC = {fextDOMred{iline,istructure}} ;
%         FD =  {fextBEAMr{iline,istructure}}  ;
%
%         if ~isempty(ELEMS_LOC)
%             %             if CUMULATIVE == 0
%             %                 Fbe_glo(ELEMS_LOC) =  {Fbe{iline,istructure}}  ;
%             %                 rRB_glo(ELEMS_LOC) = {rRB{iline,istructure}}  ;
%             %                 fextDOMred_glo(ELEMS_LOC) = {fextDOMred{iline,istructure}}  ;
%             %                 fextBEAMr_glo(ELEMS_LOC) = {fextBEAMr{iline,istructure}}  ;
%             %             else
%             Fbe_glo = AddCellsValues(Fbe_glo,ELEMS_LOC,FA,ssAA) ;
%             rRB_glo = AddCellsValues(rRB_glo,ELEMS_LOC,FB,ssBB) ;
%             fextDOMred_glo = AddCellsValues(fextDOMred_glo,ELEMS_LOC,FC,ssCC) ;
%             fextBEAMr_glo = AddCellsValues(fextBEAMr_glo,ELEMS_LOC,FD,ssDD) ;
%             % end
%         end
%
%
%
%         for iface =1:size(FORCES.TRACTIONS.LOADS,2)
%             FORCES_2_PRINT_glo(ELEMS_LOC,iface) = {FORCES_2_PRINT{iline,istructure}}  ;
%         end
%
%     end
%
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% ADDING GRAVITY
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if DATAIN.INCLUDE_GRAVITY == 1
%     for istructure =1:nstructures
%         %  elemtype = STRUCTURES(ielemtype) ;
%         ELEMS = find(MESH1D.MaterialType ==istructure); % Elements pertaining to this line and also to entity "istructure"
%
%         %  CUMULATIVE = 1 ;
%         FA = { Fgravity{istructure}}  ;
%         FB = {rRB_gravity{istructure}}  ;
%         FC = {fextDOMred_gravity{istructure}} ;
%         FD =  {fextBEAMr_gravity{istructure}}  ;
%
%         %            fextBEAMr{e,ielemtype} = fextBEAMr_traction ; %+ fextBEAMr_gravity{elemtype} ;
%         %         fextDOMred{e,ielemtype} = fextDOMred_traction ; %+ fextDOMred_gravity{elemtype} ;
%         %         Fbe{e,ielemtype} =   FtractionE ; % + Fgravity{elemtype};
%         %         rRB{e,ielemtype} =   rRB_traction ; %rRB_gravity{elemtype} ;
%
%         %if ~isempty(ELEMS)
%
%         Fbe_glo = AddCellsValues(Fbe_glo,ELEMS,FA,ssAA) ;
%         rRB_glo = AddCellsValues(rRB_glo,ELEMS,FB,ssBB) ;
%         fextDOMred_glo = AddCellsValues(fextDOMred_glo,ELEMS,FC,ssCC) ;
%         fextBEAMr_glo = AddCellsValues(fextBEAMr_glo,ELEMS,FD,ssDD) ;
%
%         %end
%
%
%     end
%
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%


%%% ASSEMBLY
% -----------------
%
%    for anod=1:nnodeE
%         a = Nod2DOF(anod,ndim) ;
%         Anod = MESH1D.CN(e,anod) ; A = Nod2DOF(Anod,ndim) ;
%         F(A) = F(A) + Fbe(a) ;
%     end

% EMPTYINDEX_1 = cellfun(@isempty,Fbe_glo) ;
% EMPTYINDEX = find(EMPTYINDEX_1==1) ;
% NEM = find(EMPTYINDEX_1==0) ;
%
% Fbe_glo(EMPTYINDEX) = {zeros(size(Fbe_glo{NEM(1)}))};
% rRB_glo(EMPTYINDEX) = {zeros(size(rRB_glo{NEM(1)}))};
% fextDOMred_glo(EMPTYINDEX) = {zeros(size(fextDOMred_glo{NEM(1)}))};
% fextBEAMr_glo(EMPTYINDEX) = {zeros(size(fextBEAMr_glo{NEM(1)}))};



end


function   fextDOMred_glo = UpdateCellForces(fextDOMred_glo,ELEMS,fextDOMred_traction)
  Facum = cell2mat(fextDOMred_glo(ELEMS)) + fextDOMred_traction ;
        [nrows ncols]= size(Facum) ;
        Facum = mat2cell(Facum,nrows,ones(1,ncols));
        fextDOMred_glo(ELEMS) = Facum ;
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
% %
% % %
% % m = nnode*ndim ; % Number of rows
% % n = 1 ;          % Number of columns
% % nzmaxLOC = size(Fbe,1)*size(Fbe,2) ;   % Maximum number of zeros (number of entries of Belems)
% % F = sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Bst
% %
% % for anod=1:nnodeE % Loop over element nodes (rows)
% %     a = Nod2DOFelem(anod,ndim,nnodeE,nelem) ;  % ROWS number (in Kelem) for node   "anod" (for all elements)
% %     %   for bnod= 1:nnodeE  % Loop over element nodes (columns)
% %     %    b = 1 ; Nod2DOF(bnod,ndim) ;
% %     Anod = MESH1D.CN(:,anod) ;  A = Nod2DOF(Anod,ndim) ;  % DOFs in the global K matrix
% %     %     Bnod = MESH1D.CN(:,bnod) ;  B = Nod2DOF(Bnod,ndim) ;
% %     %%%%%
% %     %  K(A,B) = K(A,B) + Kelem(a,b) ;
% %     %%%%%
% %     s = Fbe(a) ;
% %     s=s(:) ;
% %     % Indices "i" and "j"
% %     i = repmat(A,ndim,1);
% %     j = ones(size(i)) ;
% %     %%%%
% %
% %     F = F + sparse(i,j,s,m,n,length(s)) ;
% %     % end
% % end
% %
% % F = full(F) ;
% %
% % %
% % %    for anod=1:nnodeE
% % %         a = Nod2DOF(anod,ndim) ;
% % %         Anod = MESH1D.CN(e,anod) ; A = Nod2DOF(Anod,ndim) ;
% % %         F(A) = F(A) + Fbe(a) ;
% % %     end