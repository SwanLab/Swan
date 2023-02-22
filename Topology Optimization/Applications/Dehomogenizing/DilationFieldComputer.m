classdef DilationFieldComputer < handle
    
    properties (Access = private)
       dilation 
    end
    
    properties (Access = private)
       orientationVector
       mesh
    end
    
    methods (Access = public)
        
        function obj = DilationFieldComputer(cParams)
            obj.init(cParams)
        end
        
        function d = compute(obj)
            obj.computeDilationField();
            d = obj.dilation; 
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.orientationVector = cParams.orientationVector;
            obj.mesh  = cParams.mesh;
        end
               
        function computeDilationField(obj)
            s.fGauss = obj.computeFieldTimesDivField();
            s.mesh   = obj.mesh;
            varProb  = MinimumGradFieldWithVectorInL2(s);
            r = varProb.solve();
            s.mesh = obj.mesh;
            s.fValues = r;
            rF = P1Function(s);
            obj.dilation = rF;
        end
        
        function gradT = computeFieldTimesDivField(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CUBIC');
            a1    = obj.orientationVector{1};
            a2    = obj.orientationVector{2};            
            aDa1  = a1.computeFieldTimesDivergence(q);
            aDa2  = a2.computeFieldTimesDivergence(q);
            gradT = -aDa1.fValues - aDa2.fValues;
        end

    end

    
end