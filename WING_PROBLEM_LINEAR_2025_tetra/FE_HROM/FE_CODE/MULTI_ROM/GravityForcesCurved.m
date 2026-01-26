function [fextDOMred_glo,rRB_glo,fextBEAMr_glo,Fbe] = GravityForcesCurved(DATAIN,FORCES,ndimSP,nelem,ROTATIONS_DOMAIN_transpose,...
    DATAROM,MESH1D,nnodeE)

DATAIN = DefaultField(DATAIN,'INCLUDE_GRAVITY',1) ;
if DATAIN.INCLUDE_GRAVITY == 1
    %  if size(DATAROM{1}.GAMMAforce{1},2) == 3
    FORCES = DefaultField(FORCES,'GRAVITY', [0 -9.81  0]) ;
    
else
    FORCES.GRAVITY = [0,0,0];
end

FORCES.GRAVITY = FORCES.GRAVITY(1:ndimSP) ;
FORCES.GRAVITY = FORCES.GRAVITY(:) ;

%% ROTATION
% ---------
FG = repmat(FORCES.GRAVITY,1,nelem) ;
FG = mat2cell(FG,ndimSP,ones(1,nelem)) ;
FG =  cellfun(@mtimes,ROTATIONS_DOMAIN_transpose,FG,'UniformOutput',false) ;  % Transpose
FG = cell2mat(FG) ;

% Initialization global vectors 
fextDOMred_glo = cell(1,nelem) ;
Fbe = cell(1,nelem) ; ;
fextBEAMr_glo = cell(1,nelem)  ;
rRB_glo = cell(1,nelem)  ;
for istructure = 1:length(DATAROM)
    ELEMS = find(MESH1D.MaterialType==istructure);
    Fgravity = cell(nnodeE,1);
    for inode = 1:nnodeE
        Fgravity{inode} = DATAROM{istructure}.GAMMAforce{inode}*FG  ;
    end
    % Reduced-order force
    Fgravity= cell2mat(Fgravity) ;
    [nrows ncols]= size(Fgravity) ;
    Fgravity = mat2cell(Fgravity,nrows,ones(1,ncols));
    Fbe(ELEMS) = Fgravity ;
    %% ----------------------------------------
    fextBEAMr_gravity = DATAROM{istructure}.Lgravity*FG   ;
    [nrows ncols]= size(fextBEAMr_gravity) ;
    fextBEAMr_gravity = mat2cell(fextBEAMr_gravity,nrows,ones(1,ncols));
    fextBEAMr_glo(ELEMS) = fextBEAMr_gravity ;
    %% ----------------------------------------
    rRB_gravity= DATAROM{istructure}.BETA_gravity*FG   ;
    [nrows ncols]= size(rRB_gravity) ;
    rRB_gravity = mat2cell(rRB_gravity,nrows,ones(1,ncols));
    rRB_glo(ELEMS) = rRB_gravity ;
    %% -----------------------------------------
    fextDOMred_gravity=  DATAROM{istructure}.Dgravity*FG   ;
    [nrows ncols]= size(fextDOMred_gravity) ;
    fextDOMred_gravity = mat2cell(fextDOMred_gravity,nrows,ones(1,ncols));
    fextDOMred_glo(ELEMS) = fextDOMred_gravity ;
end