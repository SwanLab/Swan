classdef LevelSetHorizontalInclusion < LevelSetCreator
    
    properties (Access = private)
        widthH
    end

    methods (Access = public)
        
        function obj = LevelSetHorizontalInclusion(cParams)
            obj.widthH = cParams.widthH;
            obj.compute(cParams);
        end
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            m = obj.widthH;
            y0 = obj.nodeCoord;
            H = max(y0) - min(y0);
            centerY = 0.5*(max(y0) + min(y0));
            offsetY = 0.5*m*H;            
            y = y0 - centerY;
            obj.levelSet =  - (max(y./offsetY,[],2) - 1);
        end
        
    end
    
end

