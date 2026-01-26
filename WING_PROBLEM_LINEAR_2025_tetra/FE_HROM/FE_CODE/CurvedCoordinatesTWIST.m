function  COORnew  = CurvedCoordinatesTWIST(COOR,Cf1,L,angDOM)
% Transformation for twisting the slice
% -------------------------------------
if nargin == 0
    load('tmp.mat')
end
% Coordinates relative to centroid face 1 
COORrel = bsxfun(@plus,COOR,-Cf1) ;
s = COORrel(1,:) ;
angDOM_s = angDOM*s/L ;
% 
newx = COOR(1,:) ; 

newy = Cf1(2) + cos(angDOM_s).*COORrel(2,:)  - sin(angDOM_s).*COORrel(3,:) ; 
newz = Cf1(3) + sin(angDOM_s).*COORrel(2,:)  + cos(angDOM_s).*COORrel(3,:) ;

COORnew = [newx; newy; newz] ; 

 



% for idim = 1:ndim
%     if idim==1
%         COORnew(idim,:) = Cf1(idim) + cos(angDOM_s).*COORrel(1,:)  +  sin(angDOM_s).*COORrel(ind2,:) ;
%     elseif idim == indNO
%         COORnew(idim,:) = Cf1(idim) + COORrel(idim,:)  ;
%     elseif idim == ind2
%         COORnew(idim,:) = Cf1(idim) - sin(angDOM_s).*COORrel(1,:)  +  cos(angDOM_s).*COORrel(ind2,:) ;
%     end
% end


% 
% 
% ndim = size(COOR,1) ;
% 
%  
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Stretcjing along x ---
% COORnew = COOR ;
% s = COOR(1,:)-Cf1(1) ;
% angDOM_s = angDOM*s/L ;
% for idim = 1:ndim
%     if idim==1
%         COORnew(idim,:) = COOR(idim,:) +  sin(angDOM_s).*(COOR(ind2,:)-Cf1(ind2)) ;
%     elseif idim == ind2
%         COORnew(idim,:) =  Cf1(ind2) +  cos(angDOM_s).*(COOR(ind2,:)-Cf1(ind2)) ;
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

