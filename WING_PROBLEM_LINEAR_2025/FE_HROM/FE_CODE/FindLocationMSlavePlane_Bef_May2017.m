function [NslvSORT NmastSORT]= FindLocationMSlavePlane_Bef_May2017(Nmast,Nslv,COOR,TOL,INDEXplane )  

if nargin == 0
    load('tmp.mat')
end

indCOOR = setdiff(1:3,INDEXplane) ; 
COORmst = COOR(Nmast,indCOOR) ; 
COORslv = COOR(Nslv,indCOOR) ; 

%NslvSORT = zeros(size(Nslv));  % JAHOL --> We don't re-order Nslv, but
%rather redefine Nmast
NmastSORT = zeros(size(Nslv)); 

dbstop('16')
for islave=1:length(Nslv)
    coorLOC = COORslv(islave,:) ; 
    coorLOC = repmat(coorLOC,length(Nmast),1) ; 
    diffCOOR = (coorLOC-COORmst)' ; 
    normDIFF = sum(diffCOOR.*diffCOOR,1) ;
    [MINIM_c ,indiceF ]= min(normDIFF); % find(normDIFF<TOL) ; 
    if MINIM_c(1) >TOL
        %dbstop('21')
        warning('No coincidence, or multiple coincidences. Change the tolerance')
    end
    NmastSORT(islave) = Nmast(indiceF(1))  ;  
    % indiceFE is the index within Nmst ... Why not to redefine Nmast ,
    % JAHO, 3-May-2017... But how we re-order Nmast ? What happens when
    % there are several slave sets ? 
    % JAHOLaura
    %NslvSORT(islave) = Nslv(indiceF(1)) ; % What's that ??? JAHOlaura
    
end
 NslvSORT = Nslv ; 
% for imaster =1:length(Nmast)
% end