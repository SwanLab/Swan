classdef testFemComputation < handle
    
    properties (Abstract, Access = protected)
        testName
    end
    
    properties (Access = protected)
       fem
    end
    
    methods (Access = protected)
        
        function obj = testFemComputation()
           obj.computeVariableThroughFemSolver()
           obj.selectComputedVar();
        end
    end
    
    methods (Access = protected)
        
        function computeVariableThroughFemSolver(obj)
            obj.fem = FEM.create(obj.testName);
            obj.createMaterialProperties();
            obj.fem.computeVariables();
        end
        

        
    end
    
    methods (Access = private)
        
        function createMaterialProperties(obj)
            q = Quadrature.set(obj.fem.mesh.type);
            q.computeQuadrature('LINEAR');
            I = ones(obj.fem.mesh.nelem,q.ngaus);            
            p.kappa = .9107*I;
            p.mu    = .3446*I;               
            obj.fem.setMatProps(p)        
        end        
        
    end
    
    methods (Abstract, Access = protected)
        selectComputedVar(obj)
    end
    
end

