classdef SubcellsMesher_1D < SubcellsMesher_Interior
    methods (Access = protected)
        function subcells_connec = computeInteriorSubcellsConnectivities(~,~,subcell_x_value)
            subcells_connec = [1 3; 2 3];
            is_interior = all(subcell_x_value(subcells_connec) <= 0,2);
            subcells_connec = subcells_connec(is_interior,:);
        end
    end
end

