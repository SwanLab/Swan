function [MESH,ARRAY_CELLS] =  MeshCELLrepeat2Dfe(MESHrve,TRANSLATIONS,SizeArray) 
% This function returns a mesh of SizeArray(1) x SizeArray(2) of
% "generalized" elements characterized, locally, 
% by coordinates MESHrve.COOR and Connectivities CNcell 
% JAHO, 3-April-2025, UPC, CAmpus Nord
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
% --------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
 end


% At each iteration of the loop, add a new element, without worrying about possible repetitions
% It is in a second stage where such plausible repetitions are inspected. Just visit neighboring cells and removed repeated points
% Each time you remove a point, you have to update the connectivities of the tentative connectivity matrix of the current cell. 
% This implies, not only replacing the coincident nodes, but also decreasing the numbering of the remaining nodes 


% Initialization (cell 1,1)
COOR = [] ;%MESHrve.COOR ; 
nelems = prod(SizeArray) ; % Number of elements 
CN_CELL = cell(nelems,1); % Connectivities, 1 row of nodes, all 
CN = [] ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
ielem = 1; 
ARRAY_CELLS = 1:nelems ; 
ARRAY_CELLS = reshape(ARRAY_CELLS(:),SizeArray)' ; 
nelems_Acum_CONNECT = 0 ; 
for ielem_y = 1:SizeArray(2)
    for ielem_x = 1:SizeArray(1)
        % The following are the tentatives COORDINATE matrix and
        % CONNECTIVITY matrices
        NewCOOR = MESHrve.COOR + (ielem_x-1)*TRANSLATIONS{1} + (ielem_y-1)*TRANSLATIONS{2} ; 
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
            COOR =  NewCOOR   ; 
            CN  = MESHrve.CN ; 
            CN_CELL{ielem} = CN(:) ; 
          
        else
            ElemsToCheck = PATCH(IndCheck); % These are the elements whose nodes are to be checked to see if they coincide with old ones 
             NodesToCheck = unique(cell2mat(CN_CELL(ElemsToCheck))) ; % These are the existing nodes to check 
            COORtoCHECK = COOR(NodesToCheck,:) ;  % These are the coordinates of the existing nodes 
             [Idx,DIST_p] = knnsearch(COORtoCHECK,NewCOOR) ; 
             % We seek entries of DIST_p which are below a prescribed
             % tolerance (closed to zero)
             TOL = 1e-3*TRANSLATIONS{1} ; 
             IndLocREPEATED = find(DIST_p<TOL) ;  
             IndGloREPEATED = NodesToCheck(Idx(IndLocREPEATED)) ; 
             % What to do with IndLocREPEATED and IndGloREPEATED
             
             
             
             
             % What to do with IndLocREPEATED and IndGloREPEATED ?  
             % We initialize the new connectivities of the new element as
             CNnew  = MESHrve.CN(:) ; 
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