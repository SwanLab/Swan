classdef CellsClassifier < handle
    
    methods (Access = public, Static)
        
        function [full,empty,cut] = classifyCells(levelSet,connec)
            nodes = connec;
            indexes = (1:size(connec,1))';
            lsNodes = levelSet(nodes);
            phiCase = sum((sign(lsNodes)<0),2);
            
            full  = phiCase == size(nodes,2);
            empty = phiCase == 0;
            cut = indexes(~(full | empty));
        end
        
    end
    
end

