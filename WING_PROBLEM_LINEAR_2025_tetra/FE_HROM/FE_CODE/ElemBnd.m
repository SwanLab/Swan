function [CNb setBelem]= ElemBnd(CONNECTb,NODESb)
%--------------------------------------------------------------------------
% FUNCTION: ElemBnd
%
% PURPOSE:
%   Identifies and extracts boundary elements whose nodes belong entirely
%   to a given set of nodes. This is typically used to isolate the boundary
%   elements of a subdomain for integration or boundary condition imposition.
%
% USAGE:
%   [CNb, setBelem] = ElemBnd(CONNECTb, NODESb)
%
% INPUT:
%   - CONNECTb : (nElemB x nNodePerElem) connectivity matrix for boundary elements,
%                where each row corresponds to a boundary element and the columns
%                contain the global node numbers forming that element.
%   - NODESb   : List of node numbers that define the sub-boundary or subdomain.
%
% OUTPUT:
%   - CNb      : Subset of CONNECTb. Connectivity matrix of boundary elements
%                where *all* nodes belong to NODESb.
%   - setBelem : Index set of rows in CONNECTb corresponding to the selected elements.
%
% METHOD:
%   - The function loops over all columns of the boundary connectivity matrix.
%   - For each column (i.e., node in the element), it identifies which rows
%     contain nodes that are in NODESb.
%   - Using set logic and vectorization (with sorting and `unique`), it builds a
%     logical map of which elements are fully contained in the node set.
%   - The resulting `setBelem` identifies the elements whose *all* nodes are in NODESb.
%   - The alternative (legacy) loop-based implementation is retained for comparison.
%
% FEATURES:
%   - Efficient vectorized implementation using `sort`, `unique`, and `intersect`.
%   - Robust to repeated nodes in elements (e.g., quadratic elements).
%   - Includes logic to detect and handle nodes repeated more than twice
%     (see fix applied on 21-Jul-2017).
%
% APPLICATION:
%   - Frequently used in subdomain extraction, domain decomposition, patch-based
%     ROM methods, or when imposing Neumann or Dirichlet conditions on a specific face.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Version: 18-May-2017 (vectorized method)
%   Patch: 21-Jul-2017 (fix for repeated nodes in boundary elements)
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------



% CONNECTb=  Connectivity matrix for   boundary elements
% NODESb = Given set of boundary nodes
% CNb =  Connectivity matrix for  the boundary elements formed by nodes NODESb


% NEw version for removing repeated boundary elements (18-May-2017)

if nargin==0
    load('tmp2.mat')
    
%     CONNECTb = [4 7
%                 1 2
%                 2 3 
%                 2 4
%                 2 5
%                 3 6
%                 7 2
%                 ] ;
%     
%     NODESb = [2 4 7]';
    
end



METHOD =1;
NODESb = sort(NODESb) ; 

ROWS_INC_glo = zeros(size(CONNECTb)) ; 
if  METHOD ==1
    % VECTORIZED METHOD
    for i =1:size(CONNECTb,2) % Loop over nodes per element 
        % WE check which rows are candidates for being included in setBelem
       %1) Sort column i-th in ascending order --> output, the sorted vector CNsort
       % and two index vectors such that CNsort = CONNECTb(Iforward,i),  and CONNECTb(:,i) =
       % CNsort(Iback)
      [CNsort,Iforward,Iback] = sortJAHO(CONNECTb(:,i)) ;  
      % 2) Identify unique elements of CNsort. If entry a is repeated,
      % then FirstUnique(a) = LastUnique(a). 
      [CNunique FirstUnique] = unique(CNsort,'first');
      [CNunique LastUnique] = unique(CNsort,'last');     
      % In fact, FirstUnique and LastUnique serve the purpose of
      % reconstructing CNsort from CNunique. Indeed, we can say that
      % CNunique(e) occupies positions FistUnique(e):LastUnique(e) within
      % CNsort      
      % 3) Intersection between CNunique and NODESb. Output: Entries pertaining to both sets (NodesIN), and
      % an index vector OA such that CNunique(IA) = NodesIn
      [NodesIn IA] = intersect(CNunique,NODESb) ;
      % Therefore, the rows of CNunique to be included are IA. We create a
      % binary  vector with as many entries as CNsort. 
      % If CNsort(e) =1, then node "e" is to be included.  To treat
      % repeated nodes, we use the information contained in FirstUnique and
      % LastUnique
      ROWS_included = zeros(size(CNsort)) ; 
      ROWS_included(FirstUnique(IA)) = 1; 
      ROWS_included(LastUnique(IA)) = 1; 
      DiffUnique = LastUnique(IA)-FirstUnique(IA) ; % Repeated elements with more than 2 repetitions
      indRECORR = find(DiffUnique>=2) ; 
      for j = 1:length(indRECORR) 
          IAloc = IA(indRECORR(j)) ;  % Mistake detected, 21-July-2017. Before, it was IAloc = IA(j)
          ILOC  = FirstUnique(IAloc):LastUnique(IAloc) ; 
          ROWS_included(ILOC) = 1 ; 
      end
      % Go back to the original numbering
      ROWS_INC_glo(:,i) = ROWS_included(Iback) ; 
      
        
        
    end
    
    setBelem = find(prod(ROWS_INC_glo,2)==1) ;
    CNb = CONNECTb(setBelem,:) ; 
    
    
else
    setBelem = [] ;
    for ielemB=1:size(CONNECTb,1)
        %  disp(['ielemB =',num2str(ielemB),' of', num2str(size(CONNECTb,1) )])
        isSET = 1;
        for jnode=1:size(CONNECTb,2)
            jnodeG = CONNECTb(ielemB,jnode) ;
            isSET = isSET*any(jnodeG ==NODESb) ;
        end
        if isSET == 1
            setBelem = [setBelem ielemB] ;
        end
    end
    CNb = CONNECTb(setBelem,:) ;
end

end