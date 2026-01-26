function   Vnew =  FilterCandidatesInterfaceModes_AlignmentMethod(DATAIN,RotatedReactions,V,nRB,M)

if nargin == 0
    load('tmp.mat')
end
error('This function  does not work')

% Only for modes obtained by the function
% DeformModesInterface_AlignmentMethod.m

%  Each entry of V contains a set of "candidates" for being interface
%  modes. Now we have to filter out which candidates to keep.

% Diagonal matrix of interface candidates
Vcand = cell(size(V));
for i = 1:length(V)
   % if i == 1
       %   Vcand{i} = sparse(V{i}(:,nRB+1:end));  % For the first face, we only consider deformational modes  
   % else
        
    Vcand{i} = sparse(V{i});  
   % end
end
 
VcandDIAG = blkdiag(Vcand{:});
RotatedReactions = cell2mat(RotatedReactions)' ; 

% Now we define H  = Vcand*h, and search for the h  (with h^T h = ident)
% that turns H more aligned with RotatedReactions
[U,S,V]  = SVDT(VcandDIAG'*RotatedReactions) ; 
h = U*V' ; 
%%% 
H = VcandDIAG*h ; 

% 
iini = 1; 
 
for iface = 1:length(Vcand)    
  %  if iface == 1 
        % Reference face. Now we compute the intersection between H and
        % Vcand{1} and H(ROWS,:)
        VdefLOC = Vcand{iface} ;   % Candidates for def. modes
        ifin = iini +  size(VdefLOC,1) -1; 
        ROWS = [iini:ifin];  iini = ifin+1 ; 
        Hloc = H(ROWS,:) ; 
        VdefLOC = IntersectionSubspaces_M(full(Hloc),full(VdefLOC),M{iface},0,0) ; 
   % end
    end


% 