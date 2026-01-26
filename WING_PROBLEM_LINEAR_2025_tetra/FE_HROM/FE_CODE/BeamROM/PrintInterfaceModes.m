function   MSG = PrintInterfaceModes(V, DATAIN, DATA_REFMESH,NODES_faces12,MSG)

if nargin == 0
    load('tmp1.mat')
end



refFACE = 1;
COOR =DATA_REFMESH.COOR(NODES_faces12{refFACE},:) ;
CNref =  RenumberConnectivities( DATA_REFMESH.CONNECTb{refFACE},1:length(NODES_faces12{refFACE})) ;
TypeElementB = DATA_REFMESH.TypeElementB ;
posgp = [] ;
LEGENDG= ['Interface_Modes'] ;
NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),LEGENDG ];

DATAMODES = [] ;

for i = 1:size(V,2)
    V(:,i) = V(:,i)/norm(V(:,i)) ; 
end


MSG = GidPostProcessModes_dom(COOR,CNref,TypeElementB,V,posgp,NAME_MODES,DATAMODES,LEGENDG,MSG);

