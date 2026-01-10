function [MASTER SLAVES] =MasterSlavesSets(NODESpl,NODESln,NODESpnt,COOR,TOL) 

if nargin == 0
    load('tmp.mat')
end
indNODMASTER{1} = [1,1] ;  % Plane x=xmax
indNODMASTER{2} = [2,1];  % Plane y=ymax
indNODMASTER{3} = [3,1] ;  % Plane z=zmax 

indNODMASTER{4} = [1,1,2,1] ;  % Line parallel to z-axis
indNODMASTER{5} = [1,1,3,1];  % Line parallel to y-axis 
indNODMASTER{6} = [2,1,3,1] ;  % Line parallel to x-axis 

indNODMASTER{7} = 1 ;  % First corner--> This is no longer regarded as a master node ! 

%%%%% SLAVE NODES 


MASTER = cell(size(indNODMASTER)) ;
SLAVES= cell(length(indNODMASTER),7) ;

for imaster = 1:length(indNODMASTER)
    indLOC = indNODMASTER{imaster} ;
    if length(indLOC)==2
        % Plane
        MASTER{imaster} = NODESpl{indLOC(1),indLOC(2)} ;
        SLAVES{imaster,1} = NODESpl{indLOC(1),indLOC(2)+1} ;
        
        % Find exact correspondence between nodes
        % JAHOL --> Accommodate situations with more master nodes than
        % slaves nodes
        [SLAVES{imaster,1}]= FindLocationMSlavePlaneNew(MASTER{imaster},SLAVES{imaster,1},COOR,TOL,indLOC(1)) ;    
        %
        % 
        
        
    elseif length(indLOC)==4
        % Line
        MASTER{imaster} = NODESln{indLOC(1),indLOC(2),indLOC(3),indLOC(4)} ;
        SLAVES{imaster,1} = NODESln{indLOC(1),indLOC(2)+1,indLOC(3),indLOC(4)} ;
        SLAVES{imaster,2} = NODESln{indLOC(1),indLOC(2)+1,indLOC(3),indLOC(4)+1} ;
        SLAVES{imaster,3} = NODESln{indLOC(1),indLOC(2),indLOC(3),indLOC(4)+1} ;
        
        for islave = 1:3
            [ SLAVES{imaster,islave} ]= FindLocationMSlavePlaneNew(MASTER{imaster},SLAVES{imaster,islave},COOR,TOL,[indLOC(1),indLOC(3)]) ;    
        end
        
    elseif length(indLOC)==1
        % Point
        %%% All vertices are regarded as MASTER NODES !!
        MASTER{imaster} = [] ;
        restP = NODESpnt ; %setdiff(NODESpnt,NODESpnt(indLOC(1))) ; 
        for islave=1:length(restP)
           SLAVES{imaster,islave}   = restP(islave) ; 
        end
    end
end