classdef testUnfitted < test
    
    properties (Access = protected, Abstract)
        testName
        meshType
        meshIncludeBoxContour
    end
    
    properties (Access = protected)
        computation
        unfittedMesh
        oldMeshUnfitted
    end
    
    properties (Access = private)
        settings
    end
    
    methods (Access = protected)
        
    end
    
end

