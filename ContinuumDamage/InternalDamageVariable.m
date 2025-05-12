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
            tau = project(tau,obj.r.order);            
            fV = zeros(size(obj.r.fValues));
            nodesNoDamage = tau.fValues <= obj.rOld.fValues;
            fV(nodesNoDamage) = obj.rOld.fValues(nodesNoDamage);
            fV(~nodesNoDamage) = tau.fValues(~nodesNoDamage);
            obj.r.setFValues(fV);
        end

        function updateRold(obj)
            obj.rOld = obj.r;
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.rOld = project(cParams.r0,'P0');
            obj.r    = project(cParams.r0,'P0');
            obj.mesh = cParams.mesh;
        end
        
    end
    
end