function [COOR]= COOR_Quad_n9(a)
% NUMBERING CONVENTION: https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20type
% https://www.gidhome.com/documents/customizationmanual/POSTPROCESS%20DATA%20FILES/Results%20format:%20ModelName.post.res/Gauss%20Points
COOR = zeros(9,2) ;
%W = zeros(27,1) ;
% NODES 1-2-3-4: Corners  

 COOR(1,:) = [-a,-a] ; % W(1) = wLOC   ;
COOR(2,:) = [+a,-a] ; % W(2) = wLOC   ;
COOR(3,:) = [+a,+a] ; % W(3) =  wLOC  ;
COOR(4,:) = [-a,+a] ; % W(4) =  wLOC  ;
 
% NODES 5 6 7 9: Midside 
COOR(5,:) = [0,-a]  ; % W(9) = wLOC ;
COOR(6,:) = [+a,0] ; %  W(10) = wLOC ;
COOR(7,:) = [0,+a] ; % W(11) = wLOC ;
COOR(8,:) = [-a,0] ; % W(12) = wLOC ;

COOR(9,:) = [0,0] ; 

% % NODES 13-14-15-16 : Corners plane z = 0
% z = 0;
% COOR(13,:) = [-a,-a,z] ; % W(13) = wLOC ;
% COOR(14,:) = [+a,-a,z] ; % W(14) = wLOC ;
% COOR(15,:) = [+a,+a,z] ; % W(15) = wLOC ;
% COOR(16,:) = [-a,+a,z] ; % W(16) = wLOC ;
% % NODES 17-18-19-20: Midside plane z = +1
% z =  +b;
% COOR(17,:) = [0,-a,z]  ;   W(17) = wLOC ;
% COOR(18,:) = [+a,0,z] ;   W(18) = wLOC ;
% COOR(19,:) = [0,+a,z] ;   W(19) = wLOC ;
% COOR(20,:) = [-a,0,z] ;  W(20) = wLOC ;
% % NODES 21 - Center plane z = -1
% z = -b ;
% COOR(21,:) = [0,0,z] ;   W(21) = (5/9)*(8/9)^2 ;
% % NODES 22-23-24-25: Midside plane z = +1
% z =  0;    wLOC =  (5/9)*(8/9)^2  ;
% COOR(22,:) = [0,-a,z]  ;   W(22) = wLOC ;
% COOR(23,:) = [+a,0,z] ;    W(23) = wLOC ;
% COOR(24,:) = [0,+a,z] ;    W(24) = wLOC ;
% COOR(25,:) = [-a,0,z] ;   W(25) = wLOC ;
% % NODES 26 - Center plane z = +1
% z = +b ;
% COOR(26,:) = [0,0,z] ;  W(26) = wLOC ;
% % NODES 27 - Center plane z = 0
% z = 0 ;
% COOR(27,:) = [0,0,z] ;    W(27) = (8/9)^3 ;
% 
% 
% end