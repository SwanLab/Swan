classdef DilationFieldComputer < handle
    
    properties (Access = private)
       dilation 
    end
    
    properties (Access = private)
       theta 
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
        
        function plot(obj)
            figure()
            s.mesh  = obj.mesh;
            s.field = obj.dilation;
            n = NodalFieldPlotter(s);
            n.plot();
            shading interp            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.theta = cParams.theta;
            obj.mesh  = cParams.mesh;
        end
               
        function computeDilationField(obj)
            s.fGauss = obj.computeThetaGradient();
            s.mesh   = obj.mesh;
            varProb  = MinimumGradFieldWithVectorInL2(s);
            r = varProb.solve();
            obj.dilation = r;
        end
        
        function gradT = computeThetaGradient(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CUBIC');


            alpha = obj.theta;
            a1(:,1) = cos(alpha);
            a1(:,2) = sin(alpha);

            a2(:,1) = -sin(alpha);
            a2(:,2) = cos(alpha);


            s.mesh = obj.mesh;
            s.fValues = a1;            
            a1F = P1Function(s);
            a1FV = a1F.evaluate(q.posgp); 

            s.mesh = obj.mesh;
            s.fValues = a2;            
            a2F = P1Function(s);
            a2FV = a2F.evaluate(q.posgp);            

            da1 = a1F.computeDivergence(q);
            da2 = a2F.computeDivergence(q);            
            
            da1a1 = bsxfun(@times,da1.fValues,a1FV);
            da2a2 = bsxfun(@times,da2.fValues,a2FV);
          
            
            gradT = -da2a2 - da1a1;
        end

    end

    
end