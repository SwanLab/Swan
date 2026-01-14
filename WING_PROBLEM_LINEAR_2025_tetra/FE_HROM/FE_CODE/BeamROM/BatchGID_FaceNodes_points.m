function   [ TEXT]   = BatchGID_FaceNodes_points(COOR,TEXT) ;
% Write GID batch to generate a set of  points 
% See BeamROM.pdf
% JAHO, 6-Jan-2018, 4 Feb 2018
% ------------------
if nargin == 0
    load('tmp.mat')
end


%TEXT = {} ;

nnode = size(COOR,1) ;



% 
     
     for inode = 1:nnode
        TEXT{end+1} = 'Mescape Geometry Create Point' ;
        %   aaa = sort( CONNECT(inode:inode+1))   ;
        COORini = COOR(inode,:) ;
        % Check whether the nodes themselves are new or not
        TEXT{end+1} = [num2str(COORini(1)),',',num2str(COORini(2)),',',num2str(COORini(3))] ;
%         if VISIT_NODES(nnodeINI) == 1
%             TEXT{end+1} = 'old' ;
%         else
%             VISIT_NODES(nnodeINI) = 1;
%         end
        
        
     end
    
 
TEXT{end+1} = 'escape' ;



% 
% for ielem = 1:nelem
%     
%     CONNECT = [CNbF1(ielem,:) ] ;
%     for inode = 1:nnodeE
%         TEXT{end+1} = 'Mescape Geometry Create Point' ;
%         %   aaa = sort( CONNECT(inode:inode+1))   ;
%         nnodeINI =CONNECT(inode) ;
%         COORini = COOR(nnodeINI,:) ;
%         % Check whether the nodes themselves are new or not
%         TEXT{end+1} = [num2str(COORini(1)),',',num2str(COORini(2)),',',num2str(COORini(3))] ;
%         if VISIT_NODES(nnodeINI) == 1
%             TEXT{end+1} = 'old' ;
%         else
%             VISIT_NODES(nnodeINI) = 1;
%         end
%         
%         
%     end
%     
% end
% TEXT{end+1} = 'escape' ;
% 



end



