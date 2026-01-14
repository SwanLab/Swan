function  [Tnod, ElemLatSurf]= TractionForcesPrint(DATA_REFMESH,FORCES_2_PRINT,DATAIN)

if nargin == 0
    load('tmp1.mat')
end

CNlat = DATA_REFMESH.CONNECTb ;  % Connectivities boundary faces 
nfaces  = size(CNlat,2) ;   % Totla number of faces 
iface_CONNECT = 2;   % So far, only valid for slices, or two-faces joints 
iface_LATER = (iface_CONNECT+1):nfaces;    % Lateral surfaces
 
CNlat = CNlat(iface_LATER) ;  % Connectivities lateral surfaces 
                              % The boundary connectivity matrix is
                              % constructed as 
                              % CNb =  cell2mat(CNlat')
                               
FORCES_2_PRINT  =FORCES_2_PRINT(:,iface_LATER) ;  % ndomains x nfaces array 
[AA ] = cellfun(@isempty,FORCES_2_PRINT); 
ndOMS_nz =length(find(AA ==0)) ; 
Tnod = cell(ndOMS_nz,1) ; 
ElemLatSurf = Tnod ; 


ndom = size(FORCES_2_PRINT,1); 

idomNZ = 1; 
iINI = 1; 
 
for idom =1:ndom
    for iface = 1:length(CNlat)
        % What are the global indices corresponding to face "iface" of
        % domain "idom " ? 
        iFIN = iINI + size(CNlat{iface},1)-1 ; 
        if ~isempty(FORCES_2_PRINT{idom,iface})
            Tnod{idomNZ} =   FORCES_2_PRINT{idom,iface}(:) ; 
            ElemLatSurf{idomNZ} = (iINI:iFIN)' ; 
            idomNZ = idomNZ + 1; 
        end
        iINI = iFIN + 1; 
    end
end

Tnod = cell2mat(Tnod) ; 
ElemLatSurf = cell2mat(ElemLatSurf) ; 

