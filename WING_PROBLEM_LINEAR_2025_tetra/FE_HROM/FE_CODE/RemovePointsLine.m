function NODESln = RemovePointsLine(NODESln,NODESpnt) ; 

for idimA=1:size(NODESln,1)
    for jdimA = 1:2
        for idimB = 1:size(NODESln,1) 
            for jdimB = 1:2 
                LINENOD = NODESln{idimA,jdimA,idimB,jdimB} ;
                if ~isempty(LINENOD)
                    NODESln{idimA,jdimA,idimB,jdimB} = setdiff(LINENOD,NODESpnt) ; 
                end
            end
        end
    end
end