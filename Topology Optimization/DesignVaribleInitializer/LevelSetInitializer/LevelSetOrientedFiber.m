classdef LevelSetOrientedFiber < LevelSetCreator
    
    properties (Access = private)
        levelOfFibers
        fiberPosition
    end
    
    methods (Access = public)
        
        function obj = LevelSetOrientedFiber(input)
            obj.fiberPosition = input.yn;
            obj.levelOfFibers = input.levFib;
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            m = obj.levelOfFibers;
            y = obj.fiberPosition;
            period = 1/(2^m);
            phase = period/4 - mod(period/4,0.00625);
            obj.levelSet = -sin(2*pi/period*(y-phase));
        end
        
    end    

end