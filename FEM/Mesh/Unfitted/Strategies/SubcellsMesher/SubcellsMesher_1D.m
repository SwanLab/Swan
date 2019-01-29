classdef SubcellsMesher_1D < SubcellsMesher_Interior
    
    methods (Access = protected)
        
        function computeAllPossibleSubcellsInCell(obj)
            obj.subcells_connec = [1 3; 2 3];
        end
        
    end
end

