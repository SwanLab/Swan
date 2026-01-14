function [b,c,C,f,F,f_withCORNERS,F_withCORNERS]  =FaceCornersDOFS_Q4(rnodLOC,ndim,MESH,DATALOC,COMB_FACES)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/13_Q4_nocorners.mlx
%


if nargin == 0
    load('tmp.mat')
end

% STEP 1) IDENTIFY THE DISTINCT SETS OF THE DOFS  (CASE WITH CORNERS ! )
%  \item To impose periodic boundary conditions, it is necessary to divide the nodes of the interface boundary of the rectangular element in the following subsets
%
%  \begin{equation}
%   \b =  \f_1 \union  \f_2  \union \f_3  \union \f_4    \union  \c_1  \union \c_2  \union \c_3   \union \c_4
%   \end{equation}
%
%   Likewise, we denote the counterpart of the above sets in local numbering by
%
%    \begin{equation}
%   \{1,2 \ldots N\} =  \F_1 \union  \F_2  \union \F_3  \union \F_4    \union  \C_1  \union \C_2  \union \C_3   \union \C_4
%   \end{equation}
b = small2large(rnodLOC,ndim) ;  % THESE ARE THE BOUNDARY DOFS
FACES_match = DATALOC.PERIODICITY_FACES  ;

% Now we have to define c_1, c_2 ...c_4 (the corners)  and f_1, f_2, f_3 ,
% f_4 (faces without the corners )

% CORNER 1: Intersection faces
%COMB_FACES= {[1,2],[2,3],[3,4],[4,1]} ;
ndim = 2;
ncorners = 4;
%CORNERS = zeros(ncorners,1) ;
nodeCORNER = zeros(1,ncorners) ;
icorner = 1;
it_has_CORNERS = 1 ;
while icorner <= length(COMB_FACES)
    ind_FACE_A = COMB_FACES{icorner}(1) ;
    ind_FACE_B = COMB_FACES{icorner}(2) ;
    nodes_FACE_A = MESH.NODES_FACES{ind_FACE_A}  ;
    nodes_FACE_B = MESH.NODES_FACES{ind_FACE_B}  ;
    III = intersect(nodes_FACE_A,nodes_FACE_B) ;
    if ~isempty(III)
        nodeCORNER(icorner) =    III ;
    else
        it_has_CORNERS = 0  ;
        break
    end
    icorner = icorner +1 ;
end


if it_has_CORNERS == 0
    
    %\item This is the case in which the fictitious point defininig
    %the shape functions of the EIF Q4 element do not coincide with any node of the fine-scale mesh.
    ipair = 1;
    FACES_match_loc = FACES_match{ipair} ;
    nodeCORNER = cell(2,1) ;
    %  IDENTIFICATION OF THE "END NODES" OF EACH INTERFACE LINE
    for ifacesLOC = 1:length(FACES_match_loc)
        ifaces = FACES_match_loc(ifacesLOC) ;
        % CONNECTIVITIES LINE IFACES
        %---------------
        ind_CNloc =  MESH.Indexes_faces_bnd_element{ifaces};
        CNloc = MESH.CNb(ind_CNloc,:) ;
        % We are interested just in the end points of the element
        CNloc = CNloc(:,1:2);
        %
        NODESfaceLOC  = unique(CNloc(:));
        [COMMONglo,COMMONloc_1,COMMONloc_2]  =  intersect(CNloc(:,1),CNloc(:,2)) ;
        END_points = setdiff(NODESfaceLOC,COMMONglo) ;
        nodeCORNER{ifacesLOC} = END_points(:) ;
    end
    
    % REFERENCE NODE
    NODE_REF = nodeCORNER{1}(1) ;
    COORref = MESH.COOR(NODE_REF,:) ;
    % NODE other side
    distLOC = zeros(2,1) ;
    for  inode = 1:2
        nodeLOCC = nodeCORNER{2}(inode) ;
        COORloc = MESH.COOR(nodeLOCC,:) ;
        distLOC(inode) = norm(COORloc-COORref) ;
    end
    
    [dummy,indMIIN] = min(distLOC) ;
    OTHERnode =  nodeCORNER{2}(indMIIN) ;
    nodeCORNER = [NODE_REF;OTHERnode] ;
    
    
end





c = small2large(nodeCORNER,ndim) ; % Global notation d(c)
% Local notation d(b(c))
[cNEW,C,dummy2] =intersect(b,c,'stable');

c= cNEW ;

% This way  c-b(C) = 0 ;

% FACES F1,F2,F3  and F4
%
% % subSTEP 1. Remove corners
% % -------------------------
nfaces =4;
% NODES_facesNOCORN = cell(1,nfaces)
% for iface = 1:nfaces
%     nodesFACE = MESH.NODES_FACES{iface} ;
%     NODES_facesNOCORN{iface} = setdiff(nodesFACE,nodeCORNER) ;
% end

% STEPs . Sort nodes FACES such that NODE_1(i) is the opposite of NODE_3(i)



f = cell(nfaces,1) ;
F = cell(nfaces,1) ;

f_withCORNERS = cell(nfaces,1) ;
F_withCORNERS = cell(nfaces,1) ;

for imatch = 1:length(FACES_match)
    iface = FACES_match{imatch}(1) ;
    NODES_i = MESH.NODES_FACES{iface} ;
    f{iface} = small2large(NODES_i,ndim) ;
    COOR_i = MESH.PROPERTIES_FACES{iface}.COORrelA_global ;
    
    
    jface = FACES_match{imatch}(2) ;
    COOR_j = MESH.PROPERTIES_FACES{jface}.COORrelA_global ;
    NODES_j = MESH.NODES_FACES{jface} ;
    f{jface} = small2large(NODES_j,ndim) ;
    Idx_j = knnsearch(COOR_j,COOR_i) ;
    
    COOR_j_i = COOR_j(Idx_j,:) ;
    
    
    errorCOOR = norm(COOR_j_i-COOR_j,'fro') ;
    disp('---------------------------------------------')
    disp(['Checking matching interfaces = ',num2str(iface),' and ',num2str(jface)]) ;
    disp(['ERROR_match = ',num2str(errorCOOR)])
    if errorCOOR > 1e-6
        error('No periodicity in this problem')
    else
        % NODES_j(Idx_j)
        Idx_DOFs_j_i = small2large(Idx_j,ndim) ;
        f{jface} = f{jface}(Idx_DOFs_j_i) ;
    end
    
    f_withCORNERS{iface} = f{iface} ; 
    f_withCORNERS{jface} = f{jface} ;
    
    % NOW REMOVE THE CORNERS
    f{iface}  =setdiff( f{iface},c,'stable') ;
    f{jface}  =setdiff( f{jface},c,'stable') ;
    
    % Local numbering (with respect to b)
    [dummy1,F{iface},dummy2] =intersect(b,f{iface},'stable');
    [dummy1,F{jface},dummy2] =intersect(b,f{jface},'stable');
    
    
        
    % Local numbering (with respect to b)
    [dummy1,F_withCORNERS{iface},dummy2] =intersect(b,f_withCORNERS{iface},'stable');
    [dummy1,F_withCORNERS{jface},dummy2] =intersect(b,f_withCORNERS{jface},'stable');
    
    
end
