function [stressGLOv wELEM] = AverageStressOnElements(stressGLOv,wST,nelem,nstrain,ngaus)

% Rather than Gauss point stresses, we shall plot in GID average
% stresses on an element
stressGLOv = stressGLOv.*wST ;
%     strainGLOv = strainGLOv.*wST ;
stressGLOv = reshape(stressGLOv,[],nelem) ;
%    strainGLOv = reshape(strainGLOv,[],nelem) ;
wSTloc = reshape(wST,[],nelem) ;
stressELEM = zeros(nstrain,nelem);
%    strainELEM = zeros(nstrain,nelem);
wELEM = zeros(nstrain,nelem) ;
for igaus = 1:ngaus
    
    INDICES = ((igaus-1)*nstrain+1):igaus*nstrain ;
    stressELEM = stressELEM + stressGLOv(INDICES,:);
    %      strainELEM = strainELEM + strainGLOv(INDICES,:) ;
    wELEM = wELEM+wSTloc(INDICES,:) ;
end
stressELEM = stressELEM./wELEM ;
%   strainELEM = strainELEM./wELEM ;
%  strainGLOv= strainELEM(:) ;
stressGLOv = stressELEM(:) ;
wELEM = wELEM(:) ; 
%ngausLOC = 1;
