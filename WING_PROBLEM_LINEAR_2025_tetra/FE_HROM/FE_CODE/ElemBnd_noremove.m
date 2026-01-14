function [CNb setBelem]= ElemBnd_noremove(CONNECTb,NODESb)
% CONNECTb=  Connectivity matrix for   boundary elements
% NODESb = Given set of boundary nodes
% CNb =  Connectivity matrix for  the boundary elements formed by nodes NODESb

% June 22th, 2017. Version without removing repeated elements, see
% DomainDecompSD-pdf


% NEw version for removing repeated boundary elements (18-May-2017)
if nargin==0
    load('tmp.mat')
    
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

%% Remove repeated connectivities
%CONNECTb = RemoveREpeatedConnectivities(CONNECTb)  ;

METHOD =1;
NODESb = sort(NODESb) ; 

ROWS_INC_glo = zeros(size(CONNECTb)) ; 
if  METHOD ==1
    % How to do that in a vector fashion ? 
    for i =1:size(CONNECTb,2) % Loop over nodes per element 
        % WE check which rows are candidates for being included in setBelem
        
      [CNsort,Iforward,Iback] = sortJAHO(CONNECTb(:,i)) ; % Column i-th in ascending order
            [CNunique FirstUnique] = unique(CNsort,'first');
      [CNunique LastUnique] = unique(CNsort,'last');     % unique elements
      % Intersection 
      [NodesIn IA] = intersect(CNunique,NODESb) ;
      % Therefore, the rows of CNunique to b included are IA
      ROWS_included = zeros(size(CNsort)) ; 
      ROWS_included(FirstUnique(IA)) = 1; 
      ROWS_included(LastUnique(IA)) = 1; 
      DiffUnique = LastUnique(IA)-FirstUnique(IA) ; % Repeated elements with more than 2 repetitions
      indRECORR = find(DiffUnique>=2) ; 
      for j = 1:length(indRECORR) 
          IAloc = IA(j) ; 
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