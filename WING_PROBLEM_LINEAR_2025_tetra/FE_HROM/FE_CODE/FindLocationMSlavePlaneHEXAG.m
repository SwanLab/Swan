function [NslvSORT ]= FindLocationMSlavePlaneHEXAG(Nmast,Nslv,COOR,TOL,PLANES )  

% Original version: FindLocationMSlavePlane.m, see
% ALEX_THESIS_mine.pdf
%dbstop('6')
 if nargin == 0
    load('tmp.mat')
 end

 if length(Nmast) ~=length(Nslv)
   %  error('The number of master and slaves nodes is not the same')
   NslvSORT = [] ; 
   return
 end


COORmst = COOR(Nmast,:) ; 
COORslv = COOR(Nslv,:) ;  

COORmst = bsxfun(@plus,COORmst',PLANES)' ; 


NslvSORT = zeros(size(Nslv));   

for imaster=1:length(Nmast)
    coorLOC = COORmst(imaster,:) ; 
   % coorLOC = repmat(coorLOC,length(Nmast),1) ; 
    diffCOOR = bsxfun(@minus,COORslv',coorLOC') ;    
    

    normDIFF = sum(diffCOOR.*diffCOOR,1) ;
    [MINIM_c ,indiceF ]= min(normDIFF); % find(normDIFF<TOL) ; 
    if MINIM_c(1) >TOL
        %dbstop('21')
        warning('No coincidence, or multiple coincidences. Change the tolerance')
    end
%    NmastSORT(islave) = Nmast(indiceF(1))  ;  
    % indiceFE is the index within Nmst ... Why not to redefine Nmast ,
    % JAHOLaura
    NslvSORT(imaster) = Nslv(indiceF(1)) ; %
    
end
% NslvSORT = Nslv ; 
% for imaster =1:length(Nmast)
% end