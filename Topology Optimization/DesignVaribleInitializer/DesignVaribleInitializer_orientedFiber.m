classdef DesignVaribleInitializer_orientedFiber < LevelSetCreator
    
    properties (Access = private)
        levelOfFibers
        fiberPosition
    end
    
    methods
        
        function obj = DesignVaribleInitializer_orientedFiber(input)
            obj.fiberPosition = input.yn;
            obj.levelOfFibers = input.levFib;
            obj.compute(input);
            obj.computeInitialLevelSet();
        end
    end
    
    methods (Access = protected)
        
        function computeInitialLevelSet(obj)
            obj.computeLevelSet();
            obj.computeDesignVariable();
        end
        
    end    
    
    methods (Access = private)
     
        
        function computeLevelSet(obj)
            m = obj.levelOfFibers;
            y = obj.fiberPosition;
            period = 1/(2^m);
            phase = period/4 - mod(period/4,0.00625);
            obj.levelSet = -sin(2*pi/period*(y-phase));
        end
        
        function computeDesignVariable(obj)
            phi = obj.levelSet;
            switch obj.optimizerName
                case {'SLERP','HAMILTON-JACOBI'}
                    obj.x = phi;
                otherwise
                    initial_holes = ceil(max(phi,0))>0;
                    obj.x = obj.ini_design_value*ones(obj.lsSize);
                    obj.x(initial_holes) = obj.hole_value;
            end
        end
        
    end
end

