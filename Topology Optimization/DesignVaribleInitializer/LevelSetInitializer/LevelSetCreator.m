classdef LevelSetCreator < handle
    
    properties (Access = protected)
        nodeCoord
        levelSet
        lsSize
        ndim        
    end
    
    methods (Access = public)
        
        function x = getValue(obj)
            x = obj.levelSet;
        end
        
    end
    
    methods (Access = public, Static)
        
        function obj = create(d)
            factory = LevelSetFactory();
            obj     = factory.create(d);
        end
    end
    
    methods (Access = protected)
        
        function obj = compute(obj,cParams)
            obj.lsSize    = size(cParams.coord(:,1));
            obj.ndim      = cParams.ndim;
            obj.nodeCoord = cParams.coord;
            obj.computeLevelSet();
          %  obj.perturbLevelSetWhenIsZeroInNode();
        end
        
    end
    
    methods (Access = private)
       
        function perturbLevelSetWhenIsZeroInNode(obj)
           nodes = obj.levelSet == 0;
           obj.levelSet(nodes) = -1e-16;                        
        end
        
    end
    
    methods (Abstract, Access = protected)
        x = computeLevelSet(obj)
    end
    
end

