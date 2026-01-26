function   [indLINES,TEXT]   = BatchGID_FaceNodes(COOR,CNbF1,TEXT,indLINES_old) ;
% Write GID batch to generate a set of lines joining the points contained
% in CNbF1
% See BeamROM.pdf
% JAHO, 28-Dec-2017
% ------------------
if nargin == 0
    load('tmp.mat')
end
 

%TEXT = {} ;

nnode = size(COOR,1) ;
nnodeE = size(CNbF1,2) ;
nelem = size(CNbF1,1) ;
VISIT_NODES = zeros(nnode,1) ;  % Already visited nodes
INDEX_CREATED_LINES = sparse(nnode,nnode) ;

LINE_INDEX = 0;
 for ielem = 1:nelem
    SURFACES =  zeros(1,nnodeE) ; % Surfaces to be created
    
    CONNECT = [CNbF1(ielem,:) CNbF1(ielem,1)] ;
    for inode = 1:nnodeE
        TEXT{end+1} = 'Mescape Geometry Create Line' ;
        aaa = sort( CONNECT(inode:inode+1))   ;
        nnodeINI =aaa(1) ; nnodeFIN = aaa(2) ;
        % Line defined by nodeINI and  nodeFIN
        % Check whether nodeINI, nodeFIN already formed a line
        index_current_line =  INDEX_CREATED_LINES(nnodeINI,nnodeFIN) ; % == nnodeFIN) ;
        if  index_current_line==0
            % New line
            LINE_INDEX = LINE_INDEX+1;
            index_current_line = LINE_INDEX ;
            % Index new line, and nodes forming this line
            INDEX_CREATED_LINES(nnodeINI,nnodeFIN) = LINE_INDEX;
            INDEX_CREATED_LINES(nnodeFIN,nnodeINI) = LINE_INDEX;
            % Record that this line has already been created
            COORini = COOR(nnodeINI,:) ; COORfin = COOR(nnodeFIN,:);
            % Check whether the nodes themselves are new or not
            TEXT{end+1} = [num2str(COORini(1)),',',num2str(COORini(2)),',',num2str(COORini(3))] ;
            if VISIT_NODES(nnodeINI) == 1
                TEXT{end+1} = 'old' ;
            else
                VISIT_NODES(nnodeINI) = 1;
            end
            TEXT{end+1} = [num2str(COORfin(1)),',',num2str(COORfin(2)),',',num2str(COORfin(3))];
            if VISIT_NODES(nnodeFIN) == 1
                TEXT{end+1} = 'old' ;
            else
                VISIT_NODES(nnodeFIN)  =1;
            end
            
        end
        
        SURFACES(inode) = index_current_line +indLINES_old;
    end
    TEXT{end+1} = 'escape' ;
    TEXT{end+1} = 'Mescape Geometry Create NurbsSurface' ;
    TEXT{end+1} = num2str(SURFACES) ;
    TEXT{end+1} = 'escape' ;
             

    
end

indLINES= LINE_INDEX +  indLINES_old; 


