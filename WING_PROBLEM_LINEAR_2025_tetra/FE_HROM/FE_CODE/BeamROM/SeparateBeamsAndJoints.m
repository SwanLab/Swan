function [INFOLINES]= SeparateBeamsAndJoints(MESH1D)
% See BeamROM1.pdf
% OUTPUT 
% This function dissects out the geometric entities assigned to a given
% type
% INFOLINES.SUBTYPESentity  ---> Label of each subentity
% INFOLINES.NODES --> Nodes comprising each subentity (only for beams and 2-ends joints) 
% INFOLINES.ELEMENTS --> Elements comprising each subentity (only for beams and 2-ends joints) 
% INFOLINES.END_NODES --> "End" nodes of a given joint  (only for multi-connection joints)
% INFOLINES.END_ELEMENTS --> "End" elements of a given joint (only for multi-connection joints)
% INFOLINES.INTERSECT_NODES --> Intersection nodes (only for multi-connection joints)
% 
%
%
%
%
if nargin == 0
    load('tmp.mat')
end

% Types of entities (either beam or joints)
TYPESentity = unique(MESH1D.MaterialType) ;
% Loop over types of entities
nnode= size(MESH1D.COOR,1) ; % Number of 1D nodes
INFOLINES.SUBTYPESentity  = zeros(size(MESH1D.MaterialType)) ;  %
INFOLINES.NODES = cell(1,length(TYPESentity)) ;
INFOLINES.ELEMENTS = cell(1,length(TYPESentity)) ;
INFOLINES.END_NODES = cell(1,length(TYPESentity)) ;
INFOLINES.END_ELEMENTS = cell(1,length(TYPESentity)) ;
INFOLINES.INTERSECT_NODES = cell(1,length(TYPESentity)) ;

for itypeLOC = 1:length(TYPESentity)
    itype = TYPESentity(itypeLOC) ;
    % Find elements of this type
    ELEMENTS =  find(MESH1D.MaterialType == itype) ;
    % ----------------
    % Find END and INTERSECTION NODES
    % ------------------
    [BNDnodes INFOLINES.INTERSECT_NODES{itypeLOC}] = End_And_IntersectionNodes(MESH1D,ELEMENTS) ;
    
 
    
    
    if isempty(INFOLINES.INTERSECT_NODES{itypeLOC})
        % BEAMS, or JOINTS of just 2 ends
        [INDEX_SUBTYPE INFOLINES.NODES{itypeLOC} INFOLINES.ELEMENTS{itypeLOC}]= ...
            SeparateBeams(BNDnodes,MESH1D,ELEMENTS) ;
    else
        % Multiconnection JOINTS
        [INDEX_SUBTYPE,INFOLINES.END_NODES{itypeLOC}, INFOLINES.END_ELEMENTS{itypeLOC}]  = ...
            SeparateJoints(BNDnodes, INFOLINES.INTERSECT_NODES{itypeLOC},MESH1D,ELEMENTS) ;
    end
    
    INFOLINES.SUBTYPESentity(ELEMENTS)  = INDEX_SUBTYPE ;
    
end


end



function [INDEX_SUBTYPE SEQUENCE_NODES_GLO SEQUENCE_ELEMENTS_GLO]= SeparateBeams(BNDnodes,MESH1D,ELEMENTS)

CNloc = MESH1D.CN(ELEMENTS,:) ;
INDEX_SUBTYPE = zeros(size(ELEMENTS)) ;

VisitBndNodes= zeros(size(BNDnodes)) ;
isubtype = 1 ;
SEQUENCE_NODES_GLO ={} ;
SEQUENCE_ELEMENTS_GLO ={} ;

if isempty(BNDnodes)
% There may be one or several closed lines 
% How to distinghish them ?  The user has to define it as separate beams
% !!!!
error('Closed line detected. To define variations of cross-sectional area, the user has to define two different beams')
    
else
for inode = 1:length(BNDnodes)  % Loop over end nodesINTERSECT_nodes
    if VisitBndNodes(inode) == 0 % Check if the node has already been visited
        nodoINI = BNDnodes(inode) ; % Global numbering of the node
        % Look corresponding element
        [iELEM iBB] =  find(CNloc == nodoINI) ;
        % Loop until finding an "end", or a multiple connection
        SALIR = 0 ;
        SEQUENCE_NODES = [nodoINI] ;
        SEQUENCE_ELEMS = [ELEMENTS(iELEM)] ;
        while SALIR == 0
            
            INDEX_SUBTYPE(iELEM) = isubtype;
            CN_CONTIG = CNloc(iELEM,:) ;
            nodoFIN = setdiff(CN_CONTIG,nodoINI) ;
            ixxx = find(BNDnodes == nodoFIN) ;
            SEQUENCE_NODES = [SEQUENCE_NODES ; nodoFIN]  ;
            
            if ~isempty(ixxx)
                VisitBndNodes(ixxx) = 1 ;  % End of a beam/joint
                isubtype = isubtype + 1;
                break
            end
            % See which is the next element
            nodoINI = nodoFIN ;
            [iELEMSnext iBB] =  find(CNloc == nodoINI) ;
            %                 if length(iELEMSnext) >2
            %                     % This is a joint, and we are the intersection of more
            %                     % than one line
            %                     error('AMEND THIS !!!! ')
            %                     break
            %                 end
            iNEXT = setdiff(iELEMSnext,iELEM) ;
            iELEM = iNEXT ;
            SEQUENCE_ELEMS = [SEQUENCE_ELEMS ;ELEMENTS(iELEM) ] ;
            
        end
        SEQUENCE_NODES_GLO{end+1} = SEQUENCE_NODES ;
        SEQUENCE_ELEMENTS_GLO{end+1} = SEQUENCE_ELEMS ;
        
        
    end
    
    
    
end

end

end


function [INDEX_SUBTYPE  END_NODES END_ELEMENTS] = SeparateJoints(BNDnodes,INTERSECT_nodes,MESH1D,ELEMENTS) ;


CNloc = MESH1D.CN(ELEMENTS,:) ;
INDEX_SUBTYPE = zeros(size(ELEMENTS)) ;

VisitBndNodes= zeros(size(BNDnodes)) ;

END_NODES ={} ;
END_ELEMENTS = {} ; 

for inode = 1:length(INTERSECT_nodes)  % Loop over end nodesINTERSECT_nodes
    nodoINI = INTERSECT_nodes(inode) ; % Global numbering of the intersection node
    % Look corresponding elements  (at least three)
    [iELEM_int iBB] =  find(CNloc == nodoINI) ;
    
    for iBRANCHES = 1:length(iELEM_int)
        iELEM = iELEM_int(iBRANCHES) ; % First element of this branch, local numbering
        % Loop until finding an "end", or a multiple connection
        SALIR = 0 ;
        while SALIR == 0            
            INDEX_SUBTYPE(iELEM) = inode;
            CN_CONTIG = CNloc(iELEM,:) ;
            nodoFIN = setdiff(CN_CONTIG,nodoINI) ;
            ixxx = find(BNDnodes == nodoFIN) ;            
            if ~isempty(ixxx)
                END_ELEMENTS{inode}(iBRANCHES) = ELEMENTS(iELEM) ; 
                END_NODES{inode}(iBRANCHES) = nodoFIN; 
                nodoINI =  INTERSECT_nodes(inode)  ; 
                break
            end
            % See which is the next element
            nodoINI = nodoFIN ;
            [iELEMSnext iBB] =  find(CNloc == nodoINI) ;
            %                 if length(iELEMSnext) >2
            %                     % This is a joint, and we are the intersection of more
            %                     % than one line
            %                     error('AMEND THIS !!!! ')
            %                     break
            %                 end
            iNEXT = setdiff(iELEMSnext,iELEM) ;
            iELEM = iNEXT ;
            
        end
            
        
        
        
    end
    
end


end



function [BNDnodes INTERSECT_nodes] = End_And_IntersectionNodes(MESH1D,ELEMENTS)

CNloc = MESH1D.CN(ELEMENTS,:) ;
CNlocV = CNloc(:) ;
BNDnodes = [] ;
INTERSECT_nodes = [] ;
for inode = 1:length(CNlocV)
    nodeLOC  =CNlocV(inode) ;
    ix =  find(CNlocV == nodeLOC) ;
    if length(ix) ==1
        BNDnodes = [BNDnodes; CNlocV(ix)] ;
    elseif length(ix) >2
        INTERSECT_nodes = [INTERSECT_nodes ; nodeLOC] ;
    end
end
%%%%
INTERSECT_nodes = unique(INTERSECT_nodes) ;

end







