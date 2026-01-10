function BasisRrb = ConstructBasisRigidBodyORTH(COORrel)
ndim = size(COORrel,2) ; 

CENTER = sum(COORrel,1)/size(COORrel,1) ; 
COORrel = bsxfun(@plus,COORrel',-CENTER')';

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

