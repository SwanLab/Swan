function [Gi Nslave  Nmast COORsm COORcornerREF]= GmatrixSlaveMasterNEW(MASTER,SLAVES,COOR)  ;
% See explanation in document
% Plates_comp_homogenization.pdf
%dbstop('4')
if nargin == 0
    load('tmp2.mat')
end

% Number of master nodes
nmast = 0 ;
nslave = 0 ;
Nmast = [] ;
Nslave = [] ;
for i = 1:length(MASTER)
    nmast = nmast + length(MASTER{i});
    Nmast = [Nmast; MASTER{i}] ;
    for j=1:size(SLAVES,2)
        if ~isempty(SLAVES{i,j})
            nslave = nslave + length(SLAVES{i,j}) ;
            Nslave = [Nslave ;SLAVES{i,j}] ;
        end
    end
end

%diffCOOR =zeros(nslave,3) ;
COORsm.SLAVES =zeros(nslave,3) ;
COORsm.MASTER =zeros(nslave,3) ;

Gi = sparse([],[],[],nslave,nmast,nslave) ;

iniMAST = 1;
iniSLV = 1 ;

REFERENCE_CORNER = 1 ; 
icorners = 0 ; 
COORcornerREF = [] ; 
for i = 1:length(MASTER)
    finMAST = iniMAST+length(MASTER{i})-1;
    jjj = iniMAST:finMAST ;
    if ~isempty(MASTER{i})
        coorMST = COOR(MASTER{i},:) ;
    else
        % This is a corner !!!!!
        % Master coordinates --> REFERENCE CORNER 
        nodeREF = SLAVES{i,REFERENCE_CORNER} ; 
        coorMST = COOR(nodeREF,:);
        COORcornerREF = COOR(nodeREF,:); 
    end
    for j=1:size(SLAVES,2)
        
        if ~isempty(SLAVES{i,j})
            finSLV = iniSLV + length(SLAVES{i,j}) -1;
            iii = iniSLV:finSLV ;
            if ~isempty(MASTER{i})
                
                Gi = Gi+sparse(iii,jjj,1,nslave,nmast,length(iii)) ;
                
            end
            coorSLV = COOR(SLAVES{i,j},:) ;
            % diffCOOR(iniSLV:finSLV,:) = coorSLV-coorMST ;
            COORsm.SLAVES(iniSLV:finSLV,:) = coorSLV ;
            if  ~isempty(coorMST)
            COORsm.MASTER(iniSLV:finSLV,:) = coorMST ;
            end
            iniSLV = finSLV +1 ;
            
        end
        
    end
    iniMAST = finMAST + 1;
    
end


end