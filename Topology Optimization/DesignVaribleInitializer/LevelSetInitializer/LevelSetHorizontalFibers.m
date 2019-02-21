classdef LevelSetHorizontalFibers < LevelSetCreator
    
    properties (Access = private)
        levelOfFibers        
        volume
    end
    
    methods (Access = public)
        
        function obj = LevelSetHorizontalFibers(input)
            obj.levelOfFibers = input.levFib;
            obj.volume = input.volume;
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            y = obj.nodeCoord(:,2);
            m = obj.levelOfFibers;            
            nFib = 2^(m-1);
            V = obj.volume;
            yMin = min(y);
            yMax = max(y);
            L = yMax - yMin;
            T = L/(2*nFib);
            obj.levelSet = cos(2*pi*(y-yMin)/T) - cos((1-V)*pi);
        end
        
    end    

end