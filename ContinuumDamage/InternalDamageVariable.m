 classdef InternalDamageVariable < handle
    
    properties (Access = public)
        r
    end
    
    properties (Access = private)
        mesh
        rOld
    end
    
    methods (Access = public)
        
        function obj = InternalDamageVariable(cParams)
            obj.init(cParams)
        end       

        function itIs = isDamaging(obj)
            itIs = (obj.r > obj.rOld);
        end  

        function update(obj,tau)
            obj.r = max(tau,obj.rOld);
        end

        function updateRold(obj)
            obj.rOld = project(obj.r,'P0');
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.rOld = project(cParams.r0,'P0');
            obj.r    = project(cParams.r0,'P0');
            fV = 10*ones(size(obj.r.fValues));
            m = obj.r.mesh;

            % coordB = m.computeBaricenter';
            % ycoord = coordB(:,2);
            % xcoord = coordB(:,1);
            % ymiddle = (max(ycoord)-min(ycoord))/2;
            % xmiddle = (max(xcoord)-min(xcoord))/2;
            % eps = m.computeMeanCellSize/2.1;
            % isMiddleY = abs(ycoord - ymiddle) < eps;
            % isBelowMiddleX = (xcoord - (xmiddle+eps)) < 0;
            % fV(isMiddleY & isBelowMiddleX) = 100000;
            % obj.r.setFValues(fV);
            % obj.rOld.setFValues(fV);
            % obj.mesh = cParams.mesh;
        end
        
    end
    
end