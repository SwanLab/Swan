 classdef InternalDamageVariable < handle
    
    properties (Access = public)
        r
    end
    
    properties (Access = private)
        rOld
        mesh
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
            obj.r    = project(cParams.r0,'P0');
            obj.rOld = project(cParams.r0,'P0');
            obj.mesh = cParams.mesh;
        end
        
    end
    
end