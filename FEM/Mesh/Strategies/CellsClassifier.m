classdef CellsClassifier < handle
    
    methods (Access = public, Static)
        
        function [full,empty,cut] = classifyCells(levelSet,connec)
            indexes = (1:size(connec,1))';
            phiPerNodes = levelSet(connec);
            phiCase = sum((sign(phiPerNodes)<0),2);
            
            full = phiCase == size(connec,2);
            empty = phiCase == 0;
            cut = indexes(~(full | empty));
        end
        
    end
    
end

