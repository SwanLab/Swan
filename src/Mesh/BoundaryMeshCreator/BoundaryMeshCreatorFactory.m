classdef BoundaryMeshCreatorFactory < handle
    

    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.type
                case 'FromReactangularBox'
                    obj = BoundaryMeshCreatorFromRectangularBox(cParams);
                case 'FromData'
                    obj = BoundaryMeshCreatorFromData(cParams);
                case 'FromCloudPoints'
                    obj = BoundaryMeshCreatorFromCloudPoints(cParams);
            end
            
        end
        
    end
    
    
end