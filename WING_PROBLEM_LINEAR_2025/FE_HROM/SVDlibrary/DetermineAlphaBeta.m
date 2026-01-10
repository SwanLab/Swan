function [alpha beta] = DetermineAlphaBeta(A)


if iscell(A)
    %     CHECK = cellfun(@ischar,A) ;
    %     if  any(CHECK)
    %         [ aaa bbb]= find(CHECK==1) ;
    %         NAMEF = A{aaa(1),bbb(1)} ;
    %         DATA.USE_SLOW_MEMORY.ACTIVE = 1 ;
    %     end
    
    [p q]= size(A) ;
    alpha = zeros(1,p) ; beta = zeros(1,q) ;
    for iii=1:size(A,1)
        for jjj=1:size(A,2)
            
            if  ischar(A{iii,jjj})
                matObj = matfile(A{iii,jjj});
                fff = fieldnames(matObj) ;
                for kkk = 1:length(fff)
                    switch fff{kkk}
                        case 'Properties'
                        otherwise
                            var = fff{kkk} ;
                    end
                end
                info = whos(matObj,var) ;
                sizeX = info.size ;
                alpha(iii) = sizeX(1) ;
                beta(iii) = sizeX(2) ;
            else
                alpha(iii) = size(A{iii,jjj},1) ;
                beta(jjj) = size(A{iii,jjj},2) ;
            end
            
        end
    end
    
    
else
    alpha = size(A,1) ;
    beta = size(A,2) ; 
    
    
end
