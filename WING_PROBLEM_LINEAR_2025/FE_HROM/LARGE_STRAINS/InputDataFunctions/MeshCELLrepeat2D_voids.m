function [COOR,CN,ARRAY_CELLS_new] =  MeshCELLrepeat2D_voids(COORcell,TRANSLATIONS,SizeArray,VOIDS_CELLS)
% This function returns a mesh of SizeArray(1) x SizeArray(2) of
% "generalized" elements characterized, locally,
% by coordinates COORcell and Connectivities CNcell
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
ElemsAreZEro = [] ; 
ARRAY_CELLS_new = zeros(size(ARRAY_CELLS));
ielem_EFFECTIVE = 0;
for ielem_y = 1:SizeArray(2)
    for ielem_x = 1:SizeArray(1)
        % The following are the tentatives COORDINATE matrix and
        % CONNECTIVITY matrices
        if VOIDS_CELLS(ielem_y,ielem_x) == 0
            ElemsAreZEro = [ElemsAreZEro;ielem] ; 
            
        else
            ielem_EFFECTIVE = ielem_EFFECTIVE+1;
            ARRAY_CELLS_new(ielem_y,ielem_x) = ielem_EFFECTIVE ; 

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
            PATCH = setdiff(PATCH,ElemsAreZEro) ;
            
            
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
                TOL = 1e-3*norm(TRANSLATIONS{1}) ;
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
            
        end
        ielem = ielem+1;
    end
end
[ielemNEW ]= setdiff(1:nelems,ElemsAreZEro) ; 
CN = CN(ielemNEW,:) ; 

% NEW ARRAY_CELLS 
