classdef BoundaryMeshCreator < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = BoundaryMeshCreatorFactory();
            obj = f.create(cParams);
        end
        
    end
    
end