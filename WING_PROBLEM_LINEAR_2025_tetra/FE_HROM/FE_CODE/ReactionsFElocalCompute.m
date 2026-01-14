function ReactionsFINAL = ReactionsFElocalCompute(FE_NAME,DATALOC)

disp(['Reactions project = ',num2str(FE_NAME)])

if nargin == 0
    load('tmp.mat')
    
end
NameWS = [DATALOC.FOLDER,FE_NAME,'.mat'] ;

DATALOC = DefaultField(DATALOC,'ndim',3) ;

load(NameWS,'DATA_INPUT_FE','React') ;

% Nodes entities
NODES_ENTITIES = DATA_INPUT_FE.NODES_ENTITIES ;
ANG_ROTATION_TOTAL = DATA_INPUT_FE.ANG_ROTATION_TOTAL ;

ReactionsFINAL = zeros(DATALOC.ndim,length(NODES_ENTITIES))

for ientity = 1:length(NODES_ENTITIES)
    nodesLOC = NODES_ENTITIES{ientity};
    DOFS =  small2large(nodesLOC,DATALOC.ndim) ;
    ReacLOC = React(DOFS) ;
    ReacLOC = reshape(ReacLOC,DATALOC.ndim,[]) ;
    
    % Rotation
    if  ~isempty(ANG_ROTATION_TOTAL{ientity})
        ROT = ANG_ROTATION_TOTAL{ientity} ;
        ReacLOC = ROT'*ReacLOC ;
        
    end
    
    Rx = sum(ReacLOC(1,:)) ;
    Ry =  sum(ReacLOC(2,:)) ;
    Rz =  sum(ReacLOC(3,:)) ;
    disp('-----------------------------------')
    disp(['ENTITY = ',num2str(ientity)])
    disp('-----------------------------------')
    disp(['Rx = ',num2str(Rx)])
    disp(['Ry = ',num2str(Ry)])
    disp(['Rz = ',num2str(Rz)])
    
    ReactionsFINAL(1,ientity) = Rx ;
    ReactionsFINAL(2,ientity) = Ry ;
    ReactionsFINAL(3,ientity) = Rz ;
end



