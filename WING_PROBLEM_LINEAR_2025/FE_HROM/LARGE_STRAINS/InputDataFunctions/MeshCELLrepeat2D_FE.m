function [COOR,CN,ARRAY_CELLS] =  MeshCELLrepeat2D_FE(COORcell,TRANSLATIONS,SizeArray) 
% This function makes tiled copies of a given unit cell
% JAHO, 2-May-2025, Secrets by Farga, Friday, Barcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
% --------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
 %   SizeArray= [4,2] ; 
end


% At each iteration of the loop, add a new element, without worrying about possible repetitions
% It is in a second stage where such plausible repetitions are inspected. Just visit neighboring cells and removed repeated points
% Each time you remove a point, you have to update the connectivities of the tentative connectivity matrix of the current cell. 
% This implies, not only replacing the coincident nodes, but also decreasing the numbering of the remaining nodes 


% Initialization (cell 1,1)
COOR = [] ;%COORcell ; 
nelems = prod(SizeArray) ; % Number of elements 
nnodeE = size(COORcell,1) ; % Number of nodes per element
CN = zeros(nelems,nnodeE) ; 
CN_ref = 1:size(COORcell,1) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
ielem = 1; 
ARRAY_CELLS = 1:nelems ; 
ARRAY_CELLS = reshape(ARRAY_CELLS(:),SizeArray)' ; 
nelems_Acum_CONNECT = 0 ; 
for ielem_y = 1:SizeArray(2)
    for ielem_x = 1:SizeArray(1)
        % The following are the tentatives COORDINATE matrix and
        % CONNECTIVITY matrices
        NewCOOR = COORcell + (ielem_x-1)*TRANSLATIONS{1} + (ielem_y-1)*TRANSLATIONS{2} ; 
      %  NewCN = nelems_Acum_CONNECT + CN_ref ; 
        % Loop over the PATCH OF DOMAINS CENTERED AT THE STUDIED CELL 
        xPATCH = (ielem_x-1):(ielem_x+1) ; 
        iir = find(xPATCH <= SizeArray(1) &  xPATCH >0   ) ; 
        xPATCH = xPATCH(iir) ; 
        yPATCH = (ielem_y-1):(ielem_y+1) ;
        iir = find(yPATCH <= SizeArray(2) &  yPATCH >0   ) ; 
        yPATCH = yPATCH(iir) ; 
        PATCH = ARRAY_CELLS(yPATCH,xPATCH) ; 
        PATCH = PATCH(:) ; 
        % We need to only chech those cells whose label is less than ielem 
        [IndCheck ] = find(PATCH<ielem) ;
        
        if isempty(IndCheck)
            % This is the first domain 
            COOR = [COOR;NewCOOR ] ; 
            CN(ielem,:) = CN_ref ; 
          %  nelems_Acum_CONNECT = size(COOR,1) ; 
        else
            ElemsToCheck = PATCH(IndCheck); % These are the elements whose nodes are to be checked to see if they coincide with old ones 
             NodesToCheck = unique(CN(ElemsToCheck,:)) ; % These are the existing nodes to check 
            COORtoCHECK = COOR(NodesToCheck,:) ;  % These are the coordinates of the existing nodes 
             [Idx,DIST_p] = knnsearch(COORtoCHECK,NewCOOR) ; 
             % We seek entries of DIST_p which are below a prescribed
             % tolerance (closed to zero)
             TOL = 1e-3*TRANSLATIONS{1} ; 
             IndLocREPEATED = find(DIST_p<TOL) ;  
             IndGloREPEATED = NodesToCheck(Idx(IndLocREPEATED)) ; 
             % What to do with IndLocREPEATED and IndGloREPEATED ?  
             % We initialize the new connectivities of the new element as
             CNnew  = CN_ref ; 
             % Now we divide these connectibivies into two sets
             NodesToUpdate = setdiff(CNnew,IndLocREPEATED) ; 
             CNnew(NodesToUpdate)  = size(COOR,1) + (1:length(NodesToUpdate)) ; 
             CNnew(IndLocREPEATED) = IndGloREPEATED ; 
             AdditionalCOOR = NewCOOR(NodesToUpdate,:) ; 
             COOR = [COOR;AdditionalCOOR] ; 
             CN(ielem,:) = CNnew ; 
             
        end
         
        ielem = ielem+1;  
    end
end