classdef ContinuumDamageComputer < handle
    
    properties (Access = private)
        mesh
        bonudaryConditions
    end
    
    methods (Access = public)
        
        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
        end
        
        function results = compute(obj)
            K = obj.computeK();
            F = obj.computeF();
            results = obj.computeU(K,F);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end