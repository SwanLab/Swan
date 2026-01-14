function [BasisRdef_L ]= ReactionModesLocalEffectsRVE(BasisUdef_L,BasisRdef_L,f,TOL_SINGULAR_VALUES_Hqr)

% Copy of ReactionModesLocalEffectsNEW
if nargin == 0
    load('tmp1.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = BasisUdef_L(f,:) ;
F = BasisRdef_L(f,:) ;
%
H  = U'*F ; % Work done by each reaction force F(:,1), F(:,2) ...
% Let us calculate it in a dimensionless fashion (dividing it by the product of their norm)
normH = sqrt(sum(H.*H,1)) ;
normF = sqrt(sum(F.*F,1)) ;
normU = norm(U,'fro');
% Therefore
WORK = normH./(normF*normU)  ;

% TOLERANCE for including
TOL = 0.1;
IND_FORCES = find(abs(WORK)>TOL) ;
if isempty(IND_FORCES)
    % NO REACTION FORCES are INCORPORATED
    BasisRdef_L = [] ;
    IMODES_INCLUDE = [] ; 
else   
     % Now we have to decide wich columns
     % To this end, we first sort the columns of "" in descending order
     [~,INXmodes] = sort(WORK,'descend')  ;
     INXmodes = 1:length(WORK) ; 
    % Determining number of reaction modes
    %b -----------------------------------
    imodeLOC = 1;
    
    MODES_INCLUDE =[] ;
    IMODES_INCLUDE = [] ;
    while  imodeLOC <=size(F,2)
        imode = INXmodes(imodeLOC) ; 
        NEW_MODES = [MODES_INCLUDE,F(:,imode) ] ;
        HqrT = NEW_MODES'*U;
        SSVAL = svd(HqrT) ;
        if imodeLOC == 1
            ratioSV =1 ; % SSVAL(1) ;
        else
            ratioSV = SSVAL(end)/SSVAL(end-1) ;
        end
        if ratioSV >= TOL_SINGULAR_VALUES_Hqr
            MODES_INCLUDE = NEW_MODES ;
            IMODES_INCLUDE(end+1) = imode ;
        end
        imodeLOC = imodeLOC + 1;
        
        
    end
    
     BasisRdef_L = BasisRdef_L(:,IMODES_INCLUDE);
    
    
    
end


