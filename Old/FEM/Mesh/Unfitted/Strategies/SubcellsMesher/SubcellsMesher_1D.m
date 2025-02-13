classdef SubcellsMesher_1D < SubcellsMesher_Interior_Abstract
    
    methods (Access = public)
        
        function obj = SubcellsMesher_1D(cParams)
           obj.init(cParams); 
        end
        
    end    
    
    methods (Access = protected)
        
        function computeAllPossibleSubcellsInCell(obj)
            obj.subcells_connec = [1 3; 2 3];
        end
        
    end
end

