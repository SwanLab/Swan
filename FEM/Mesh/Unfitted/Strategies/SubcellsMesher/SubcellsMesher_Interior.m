classdef SubcellsMesher_Interior < SubcellsMesher_Interior_Abstract
    
  
    methods (Access = public)
        
        function obj = SubcellsMesher_Interior(cParams)
           obj.init(cParams); 
        end
        
    end
    
    methods (Access = protected)
        
        function computeAllPossibleSubcellsInCell(obj)
            if size(obj.coord_iso,2) == 1
                [~,I] = sort(obj.coord_iso);
                connec = [I circshift(I,-1)];
                connec(end,:) = [];
                obj.subcells_connec = connec;
            else
                obj.subcells_connec = obj.computeDelaunay(obj.coord_iso);
            end
        end
        
    end
    
end

