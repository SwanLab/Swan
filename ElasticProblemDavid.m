classdef ElasticProblemDavid < handle
    
    properties (Access = public)
        uFun
        sigma
    end
    
    properties (Access = private)
            





    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ElasticProblemDavid()
            obj.init()
            obj.createMesh()
            obj.createMaterial()
            obj.createBoundaryConditions()
            obj.displacementFunction()
            obj.createStiffness()
            obj.createForce()
            obj.computeNewDisplacement()
            obj.computeStress()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
    end
    
end