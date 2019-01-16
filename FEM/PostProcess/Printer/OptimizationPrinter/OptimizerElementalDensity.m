classdef OptimizerElementalDensity < handle
    
    properties (Access = protected, Abstract)
       fields 
    end
    
    methods (Access = protected)
        
        function createTopOptFields(obj,x,cost,constraint)
            ls = x;
            phyPr = cost.ShapeFuncs{1}.getPhysicalProblem();
            filter = FilterP0(ls,phyPr);
            dens = filter.getDens0();
            obj.fields.dens = dens;
        end
    end
    
    
end