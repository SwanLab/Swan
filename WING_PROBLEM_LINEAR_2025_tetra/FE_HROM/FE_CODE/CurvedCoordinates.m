function  COORnew  = CurvedCoordinates(COOR,Cf1,L,angDOM,ISBEAM)
% Transformation
% --------------
% Rotation matrix face 3  (relative to reference system of face 1)
% R3 = [cos(angDOM_t), 0 , sin(angDOM_t)
%      0            1        0
%      -sin(angDOM_t)  0     cos(angDOM_t)] ;
if nargin == 4
    ISBEAM =0 ; 
end
ndim = size(COOR,1) ;


if ndim == 2 | ISBEAM ==1 ; 
    ind2 = 2;
    indNO = 3;
else
    ind2 = 3;
    indNO = 2 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stretcjing along x ---
COORnew = COOR ;
%COORnew(1,:) = COOR(1,:) +  angDOM*(COOR(ind2,:)-Cf1(ind2)).*(COOR(1,:)-Cf1(1) )/L ;
s = COOR(1,:)-Cf1(1) ;
angDOM_s = angDOM*s/L ;
for idim = 1:ndim
    if idim==1
        COORnew(idim,:) = COOR(idim,:) +  sin(angDOM_s).*(COOR(ind2,:)-Cf1(ind2)) ;
    elseif idim == ind2
        COORnew(idim,:) =  Cf1(ind2) +  cos(angDOM_s).*(COOR(ind2,:)-Cf1(ind2)) ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COORrelORIG = bsxfun(@plus,COOR,-Cf1(:)) ;
COORrel = bsxfun(@plus,COORnew,-Cf1(:)) ;
COORnew = zeros(size(COORrel)) ;
s = COORrelORIG(1,:) ;
angDOM_s = angDOM*s/L ;
for idim = 1:ndim
    if idim==1
        COORnew(idim,:) = Cf1(idim) + cos(angDOM_s).*COORrel(1,:)  +  sin(angDOM_s).*COORrel(ind2,:) ;
    elseif idim == indNO
        COORnew(idim,:) = Cf1(idim) + COORrel(idim,:)  ;
    elseif idim == ind2
        COORnew(idim,:) = Cf1(idim) - sin(angDOM_s).*COORrel(1,:)  +  cos(angDOM_s).*COORrel(ind2,:) ;
    end
end
