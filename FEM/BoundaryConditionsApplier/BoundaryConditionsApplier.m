classdef BoundaryConditionsApplier < handle

    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = BoundaryConditionsApplierFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public, Abstract)
        fullToReducedMatrix(obj) 
        fullToReducedVector(obj)
        reducedToFullVector(obj)
    end
   
end

