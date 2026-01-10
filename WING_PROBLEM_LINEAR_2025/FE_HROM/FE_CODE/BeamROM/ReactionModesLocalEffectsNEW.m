function [BasisRdef_L_new,TEXT ]= ReactionModesLocalEffectsNEW...
    (BasisUdef_L,BasisRdef_L,f,TOL_SINGULAR_VALUES_Hqr,SinvVal_Rdef,nrigidB,DATAIN)

if nargin == 0
    load('tmp2.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = BasisUdef_L(f,:) ;
F= BasisRdef_L(f,:) ;
SinvVal_Rdef = SinvVal_Rdef/SinvVal_Rdef(1) ;
Fclas = bsxfun(@times,BasisRdef_L(f,:)',SinvVal_Rdef)' ;
%
H  = U'*Fclas ; % Work done by each reaction force F(:,1), F(:,2) ...
% Let us calculate it in a dimensionless fashion (dividing it by the product of their norm)
normH = sqrt(sum(H.*H,1)) ;
%normF = sqrt(sum(F.*F,1)) ;

% Therefore
WORK = normH ;

% TOLERANCE for including
TEXT = {} ; 
TOL = 1e-6 ;
IND_FORCES = find(abs(WORK)>TOL) ;

% Set it by default in zero. Numerical experience shows that this far
% better. 
DATAIN = DefaultField(DATAIN,'SORT_REACTION_MODES_ACCORDING_TO_WORK_DONE',0) ; 

if isempty(IND_FORCES)
    % NO REACTION FORCES are INCORPORATED
    BasisRdef_L = [] ;
    IMODES_INCLUDE = [] ;
else
    % Now we have to decide wich columns
    % To this end, we first sort the columns of "" in descending order
    if DATAIN.SORT_REACTION_MODES_ACCORDING_TO_WORK_DONE == 1
        [~,INXmodes] = sort(WORK,'descend')  ;
    else
        INXmodes = 1:length(WORK) ;
    end
    %   ;
    % Determining number of reaction modes
    %b -----------------------------------
    imodeLOC = 1;
    
    MODES_INCLUDE =[] ;
    IMODES_INCLUDE = [] ;
    REJECTED = [] ;
    ratioSV_REJECTED = [] ;
    while  imodeLOC <=size(F,2) && length(IMODES_INCLUDE)< size(BasisUdef_L,2)
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
        else
            REJECTED =[REJECTED, imode] ;
            ratioSV_REJECTED = [ratioSV_REJECTED,ratioSV] ;
            
            TEXT{end+1} = ['React. Mode ',num2str(imode),' rejected'] ; 
                        TEXT{end+1} = ['ratioSV = ',num2str(ratioSV)] ; 

            
        end
        imodeLOC = imodeLOC + 1;
        
        
    end
    BasisRdef_L_new = BasisRdef_L(:,IMODES_INCLUDE);
    
    if size(BasisRdef_L_new,2) < nrigidB
        % We have to complete the basis matrix
        % -----------------------------------
        SinvVal_Rdef = SinvVal_Rdef(REJECTED) ;
        ratioSV_REJECTED = ratioSV_REJECTED.*SinvVal_Rdef' ; 
        warning('Not complete basis matrix for reactions... Augmenting it...')
        [~,IX ]= sort(ratioSV_REJECTED,'descend') ;
        IX = IX(1:(nrigidB- size(BasisRdef_L_new,2))) ;
        imode = INXmodes(REJECTED(IX)) ;
        IMODES_INCLUDE = [IMODES_INCLUDE, imode] ;
        BasisRdef_L_new = BasisRdef_L(:,IMODES_INCLUDE);
    end
    
    
end

TEXT{end+1} = ['Included reaction modes'] ; 
TEXT{end+1} = [num2str(sort(IMODES_INCLUDE))] ; 


TEXT{end+1} = ['**************************************************+'];


