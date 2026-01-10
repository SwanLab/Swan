function [CNnew,TRANSF_COORD,Vall_rot,PERMUT_chosen] = ParentDomainSearchEIFEsupp(EIFEoper_all,COOR,CNold,TypeElement,DATA,nnodeE_support)
% -------------------------------------------------------------------------
%  Given a EIF element formed by nodes CNold, and the properties of
%  a set of parent domain candidates EIFEoper_all, this function searches
%  for the parent domain that minimizes the dilatational component of the
%  transformation of coordinates from the physical domain to the parent
%  domain
% OUTPUT:
% CNnew: New order of the nodes of the EIF element (permutation)
% TRANSF_COORD =
%
%   struct with fields:
%
%        ROTATION: [2×2 double]
%     TRANSLATION: [2×1 double]
%     SCALEFACTOR: 0.0500
%           detJe: 0.0025
%
% Vall_rot --> Rotated interface displacement modes
% ----------------------------------------------------------
%
% JAHO, 10-March-2023/11-March-2023
% -------------------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
end
% ---------------------------------------------------------------------
% DETERMINING PARENT DOMAIN FROM A SET OF POTENTIAL CANDIDATES
% ---------------------------------------------------------------------
icand = 1;  % EIFEoper_all(1), EIFEoper_all(2) ....
% When data is provided this way, we have to check which transformation
% yields zero dilatation
DATA = DefaultField(DATA,'ToleranceDeformationalPartTransformationParentDomain',1e-5) ;
while icand <=length(EIFEoper_all)
    
    
    % PERMUTATION CONNECTIVITIES
    
    TRANSF_COORD_perm = cell(1,length(DATA.PERMUT));
    normDmin = 1e20 ; 
    for   iPERM = 1:length(DATA.PERMUT)
        %   disp(['PERM = ',num2str(iPERM)])
        CNnew = CNold(DATA.PERMUT{iPERM}) ;
        CNnew_input = CNnew(1:nnodeE_support) ; % 26-Apr-2204
        Xe = COOR(CNnew_input,:)' ;
        [TRANSF_COORD_perm{iPERM},normD] = PolarDecompEIFEperm(Xe,EIFEoper_all(icand),DATA) ;
        normDmin = min(normDmin,normD) ; 
    end
    
    [TRANSF_COORD,PERMUT_chosen,CNnew] = CriterionChooseParent(TRANSF_COORD_perm,CNold,DATA,DATA.PERMUT) ;
    
    
  
    
    if  isempty(TRANSF_COORD)
        icand = icand + 1;
    else
        %EIFEoper = EIFEoper_all(icand);
        TRANSF_COORD.IndexParentDomain = icand ;
        [ndim,~] = size(Xe) ;
        Vall_rot = zeros(size( EIFEoper_all(icand).MODES.Vall)) ;
        
        %  UNIFORM_SCALING_REFERENCE_ELEMENT= 1;
        
        %         if  DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1
        %             FactorV = (1/TRANSF_COORD.SCALEFACTOR) ;
        %             % THIS IS TO REMIND US THAT THIS SCALING FACTOR ACTUALLY
        %             % APPEARS IN DERIVING THE FE B-matrix. We incorporate it here
        %             % because this version (11-March-2023) can only accurately handle uniform
        %             % dilations
        %         else
        %             FactorV = 1;
        %         end
        
        for imodes = 1:size( EIFEoper_all(icand).MODES.Vall,2)
            % Rotation of the interface modes
            % LOCV = FactorV*(TRANSF_COORD.ROTATION*reshape(EIFEoper_all(icand).MODES.Vall(:,imodes),ndim,[])) ;
            % disp('Temporal...Erase it')
            
            
            LOCV = (TRANSF_COORD.ROTATION'*reshape(EIFEoper_all(icand).MODES.Vall(:,imodes),ndim,[])) ;
            Vall_rot(:,imodes) = LOCV(:) ;
        end
        
        
        
        break
    end
    
    
    
    
end
if icand >length(EIFEoper_all)
  
    error(['THIS ELEMENT HAS NO PARENT DOMAIN; maximum gap =  ',num2str(normDmin),  '(you may fix it by setting  DATAOUT.ToleranceDeformationalPartTransformationParentDomain above this value   )'])
end