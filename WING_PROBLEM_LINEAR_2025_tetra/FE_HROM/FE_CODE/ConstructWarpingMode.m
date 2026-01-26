function BasisRrb = ConstructWarpingMode(COORrel)

 

ndim = size(COORrel,2) ;

if ndim == 2
  % error('Option not valid')
   BasisRrb = [] ; 
else
    
    BasisRrb =  zeros(prod(size(COORrel)),1) ;
%     BasisRrb(1:3:end) = COORrel(:,3).*COORrel(:,2) ;
%     BasisRrb(2:3:end) = 0;
%     BasisRrb(3:3:end) = 0;
    
    BasisRrb(1:3:end) = -COORrel(:,3).*COORrel(:,2);
    BasisRrb(3:3:end) = +COORrel(:,1).*COORrel(:,2);
    
    
%     % Rotation modes
%     BasisRrb(2:3:end,4) = +COORrel(:,3);
%     BasisRrb(3:3:end,4) = -COORrel(:,2);
%     
%     BasisRrb(1:3:end,5) = -COORrel(:,3);
%     BasisRrb(3:3:end,5) = +COORrel(:,1);
%     
%     BasisRrb(1:3:end,6) = -COORrel(:,2);
%     BasisRrb(2:3:end,6) = +COORrel(:,1);
    
end
% 
% DATAIN = DefaultField(DATAIN,'ORTHOGONAL_RIGID_BODY_MODES',0) ;
% 
% if DATAIN.ORTHOGONAL_RIGID_BODY_MODES== 1
%     BasisRrb = orth(BasisRrb)  ;
% elseif DATAIN.ORTHOGONAL_RIGID_BODY_MODES== 2
%     NNN = sqrt(sum(BasisRrb.^2,1)) ;
%     for imode = 1:size(BasisRrb,2)
%         BasisRrb(:,imode) = BasisRrb(:,imode)/NNN(imode) ;
%     end
% end