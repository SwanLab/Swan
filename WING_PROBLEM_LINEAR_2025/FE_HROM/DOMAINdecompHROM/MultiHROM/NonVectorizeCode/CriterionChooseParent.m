function [TRANSF_COORD,PERMUT_chosen,CNnew] = CriterionChooseParent(TRANSF_COORD_perm,CNold,DATA,PERMUT)

if nargin == 0
    load('tmp.mat')
end
TRANSF_COORD =[] ; PERMUT_chosen = [] ; CNnew =[] ;
%   TRANSF_COORD_perm = cell(1,length(PERMUT));
%     Q11 = zeros(length(PERMUT),1) ;
%     for   iPERM = 1:length(PERMUT)
%         CNnew = CNold(PERMUT{iPERM}) ;
%         Xe = COOR(CNnew,:)' ;
%         [TRANSF_COORD_perm{iPERM}] = PolarDecompEIFEperm(Xe,EIFEoper_all(icand)) ;
%         if ~isempty(TRANSF_COORD_perm{iPERM})
%             Q11(iPERM) = TRANSF_COORD_perm{iPERM}.ROTATION(1,1) ;
%         end
%     end

Q11 = zeros(length(TRANSF_COORD_perm),1) ;
Q33= zeros(length(TRANSF_COORD_perm),1) ;

for   iPERM = 1:length(TRANSF_COORD_perm)
    
    if ~isempty(TRANSF_COORD_perm{iPERM})
        if  size( TRANSF_COORD_perm{iPERM}.ROTATION,2) == 2
            Q11(iPERM) = atand(TRANSF_COORD_perm{iPERM}.ROTATION(1,2)/TRANSF_COORD_perm{iPERM}.ROTATION(1,1)) ;
            
            ang_cos = acosd(TRANSF_COORD_perm{iPERM}.ROTATION(1,1)) ;
            ang_sin = asind(TRANSF_COORD_perm{iPERM}.ROTATION(1,2)) ;
            
            Q11(iPERM) = real(min(ang_sin,ang_cos)) ;
            
            
        else
            Q11(iPERM) = TRANSF_COORD_perm{iPERM}.ROTATION(1,1) ;
            Q33(iPERM) = TRANSF_COORD_perm{iPERM}.ROTATION(3,3) ;
        end
    end
end


switch DATA.CriterionChooseParentDomain
    case 'MAX_Q_11'
        IndNonEmpty = find(~cellfun(@isempty, TRANSF_COORD_perm)) ;
        if  ~isempty(IndNonEmpty)
            [~,index_choose] = max(Q11(IndNonEmpty)) ;
            index_choose = IndNonEmpty(index_choose) ;
            TRANSF_COORD = TRANSF_COORD_perm{index_choose(1)} ;
            PERMUT_chosen = PERMUT{index_choose} ;
            CNnew = CNold(PERMUT_chosen) ;
        else
            TRANSF_COORD = [] ;
        end
        
    case 'MIN_Q_11'
        IndNonEmpty = find(~cellfun(@isempty, TRANSF_COORD_perm)) ;
        if  ~isempty(IndNonEmpty)
            [~,index_choose] = min(Q11(IndNonEmpty)) ;
            index_choose = IndNonEmpty(index_choose) ;
            TRANSF_COORD = TRANSF_COORD_perm{index_choose(1)} ;
            PERMUT_chosen = PERMUT{index_choose} ;
            CNnew = CNold(PERMUT_chosen) ;
        else
            TRANSF_COORD = [] ;
        end
        
    case 'MAX_Q_33'
        IndNonEmpty = find(~cellfun(@isempty, TRANSF_COORD_perm)) ;
        if  ~isempty(IndNonEmpty)
            [~,index_choose] = max(Q33(IndNonEmpty)) ;
            index_choose = IndNonEmpty(index_choose) ;
            TRANSF_COORD = TRANSF_COORD_perm{index_choose(1)} ;
            PERMUT_chosen = PERMUT{index_choose} ;
            CNnew = CNold(PERMUT_chosen) ;
        else
            TRANSF_COORD = [] ;
        end
        
    otherwise
        error('Option not implemented')
        
end