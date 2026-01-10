function [Gi Nslave  Nmast diffCOOR]= GmatrixSlaveMaster(MASTER,SLAVES,COOR)  ;
% See FEM_theory.pdf, page 86
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

ndim = size(COOR,2) ; 

diffCOOR =zeros(nslave,ndim) ;

Gi = sparse([],[],[],nslave,nmast,nslave) ;

iniMAST = 1;
iniSLV = 1 ;
for i = 1:length(MASTER)
    finMAST = iniMAST+length(MASTER{i})-1;
    jjj = iniMAST:finMAST ;
    if ~isempty(MASTER{i})
        coorMST = COOR(MASTER{i},:) ;
    else
        coorMST = 0 ;
    end
    for j=1:size(SLAVES,2)
        
        if ~isempty(SLAVES{i,j})
            finSLV = iniSLV + length(SLAVES{i,j}) -1;
            iii = iniSLV:finSLV ;
            if ~isempty(MASTER{i})
                
                Gi = Gi+sparse(iii,jjj,1,nslave,nmast,length(iii)) ;
                
            end
            coorSLV = COOR(SLAVES{i,j},:) ;
            diffCOOR(iniSLV:finSLV,:) = coorSLV-coorMST ;
            
            iniSLV = finSLV +1 ;
        end
        
    end
    iniMAST = finMAST + 1;
    
end


end