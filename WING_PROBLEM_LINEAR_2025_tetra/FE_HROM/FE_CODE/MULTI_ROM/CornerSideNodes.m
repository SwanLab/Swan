function [NODES_CORNERS,NODES_SIDES,NODES_FACES] = CornerSideNodes(DATA3D)


%%% CONTINUOUS STRUCTURES ---
% ---------------------------
%  Check whether there is intersection between DOFs of faces
%  4-1, 1-2, 2-3, 3-4
%  ***  ***  ***  ***
%  C1   C2   C3   C4
% -----------------------------------------------------------
% OUTPUT
 NODES_CORNERS = [] ;
 NODES_SIDES = [] ;
 NODES_FACES = [] ; 
% -------------------------
PAIRS = {[4,1],[1,2],[2,3],[3,4]} ;

% CORNER NODES
% --------------
for ipair = 1:length(PAIRS)
    iFACE1 = PAIRS{ipair}(1) ;
    iFACE2 = PAIRS{ipair}(2) ;
    [CORNER  IXX]=  intersect(DATA3D.NODES_FACES{iFACE1},DATA3D.NODES_FACES{iFACE2}) ;
    if ~isempty(CORNER)
       NODES_CORNERS{ipair} = CORNER ;
    end
    %%
    
end

% SIDE NODES
% -----------
if  ~isempty(NODES_CORNERS)
    PAIRScorners = {[1,2],[2,3],[3,4],[4,1]} ;
    for iside = 1:length(PAIRScorners)
        iCORNER1 = PAIRScorners{iside}(1) ;
        iCORNER2 = PAIRScorners{iside}(2) ;
        CORNERS = [NODES_CORNERS{iCORNER1}',NODES_CORNERS{iCORNER2}'] ;
        SIDE = setdiff(DATA3D.NODES_FACES{iside},CORNERS) ;
        NODES_SIDES{iside} = SIDE ;
        %%% Now we re-arrange the entries of DATA3D.NODES_FACES{iface} as follows 
       % ---------------------------------------------------------------------------
       %  FACE 1: CORNER1 - CORNER2- SIDE 1 
       % ------  
       if iside == 1 | iside == 2 
       DATA3D.NODES_FACES{iside} = [CORNERS,NODES_SIDES{iside}']';   
       else
        iCORNER1 = PAIRScorners{iside}(2) ;
        iCORNER2 = PAIRScorners{iside}(1) ;
        CORNERS = [NODES_CORNERS{iCORNER1}',NODES_CORNERS{iCORNER2}'] ;
         DATA3D.NODES_FACES{iside} = [CORNERS,NODES_SIDES{iside}']';   
       end
 
    end


 

end

NODES_FACES =   DATA3D.NODES_FACES ; 

