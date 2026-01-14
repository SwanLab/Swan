function DOFSeliminate = DOFScolumn(ndim,MESH2D)

[ndimMAX faceMAX]= max(ndim) ; 
ndimX = ndim(1) ; 
ndimY = ndim(2) ; 
% Which are the faces with less modes ?  --> faceMIN 
% -------------------------------------   -----------
if faceMAX == 1
    faceMIN = [2 4] ;
else
    faceMIN = [1 3] ; 
end
% -------------------------------------------------
% DOFs to eliminate 
% -----------------------------------------
NODESmin = unique(MESH2D.CN(:,faceMIN)) ; 
DOFSmin = small2large(NODESmin,ndimMAX) ; 
DOFSmin = reshape(DOFSmin,ndimMAX,[]) ; 
idofADD = ndim(faceMIN(1))+1:ndimMAX ; 
DOFSeliminate = DOFSmin(idofADD,:) ; 