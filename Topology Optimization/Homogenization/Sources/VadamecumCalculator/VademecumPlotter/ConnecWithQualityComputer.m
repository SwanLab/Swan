classdef ConnecWithQualityComputer < handle
    
    properties (Access = private)
        x
        y
        badConnec
    end
    
    methods (Access = public)
        
        function obj = ConnecWithQualityComputer(cParams)
            obj.init(cParams)            
        end
        
        function connec = compute(obj)
            isQ = obj.isNotBadQuality | obj.isNotLargeDensityElement();    
            connec = obj.badConnec(isQ,:);     
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.x      = cParams.x;
            obj.y      = cParams.y;
            obj.badConnec = cParams.connec;
        end
        
        function m = createMesh(obj)
            s.coord = [obj.x,obj.y];
            s.connec = obj.badConnec;
            m = Mesh.create(s);            
        end
        
        function itIsNot = isNotBadQuality(obj)
            m = obj.createMesh();
            qua = m.computeElementQuality';            
            itIsNot = ~(qua < 0.1);            
        end
        
        function isNot = isNotLargeDensityElement(obj)
            c = obj.badConnec;
            yV = obj.y;
            yMean = 1/3*(yV(c(:,1)) + yV(c(:,2)) +yV(c(:,3)));            
            isNot = (yMean > 0.8);
        end

    end
    
end