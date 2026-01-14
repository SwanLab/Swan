function BasisRrb = ConstructBasisRigidBody_centroid(COORabs,DATAIN)

if nargin == 1
    DATAIN.ORTHOGONAL_RIGID_BODY_MODES = 0 ;
end

nnode = size(COORabs,1) ;
Centoid = sum(COORabs,1)/nnode; % Centroid
COORrel = bsxfun(@minus,COORabs',Centoid')'; % Coordinates relative to centroid
ndim = size(COORabs,2) ; 
if ndim == 2
    BasisRrb =  zeros(prod(size(COORrel)),3) ;
    BasisRrb(1:ndim:end,1) = 1;
    BasisRrb(2:ndim:end,2) = 1;
    BasisRrb(1:ndim:end,3) = COORrel(:,2);
    BasisRrb(2:ndim:end,3) = -COORrel(:,1);  % or the other way around ?
else
    
    BasisRrb =  zeros(prod(size(COORrel)),6) ;
    BasisRrb(1:3:end,1) = 1 ;
    BasisRrb(2:3:end,2) = 1;
    BasisRrb(3:3:end,3) = 1;
    % Rotation modes
    BasisRrb(2:3:end,4) = COORrel(:,3);
    BasisRrb(3:3:end,4) = -COORrel(:,2);
    
    BasisRrb(1:3:end,5) = -COORrel(:,3);
    BasisRrb(3:3:end,5) = COORrel(:,1);
    
    BasisRrb(1:3:end,6) = COORrel(:,2);
    BasisRrb(2:3:end,6) = -COORrel(:,1);
    
end

DATAIN = DefaultField(DATAIN,'ORTHOGONAL_RIGID_BODY_MODES',0) ;

if DATAIN.ORTHOGONAL_RIGID_BODY_MODES== 1
    BasisRrb = orth(BasisRrb)  ;
elseif DATAIN.ORTHOGONAL_RIGID_BODY_MODES== 2
    NNN = sqrt(sum(BasisRrb.^2,1)) ;
    for imode = 1:size(BasisRrb,2)
        BasisRrb(:,imode) = BasisRrb(:,imode)/NNN(imode) ;
    end
end