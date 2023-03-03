classdef ScalarProduct < handle
    
    properties (Access = public)
        epsilon
        Ksmooth
        Msmooth
    end
    
    properties (Access = private)
       nVariables
       mesh
       field
    end
    
    methods (Access = public)
        
        function obj = ScalarProduct(cParams)
            obj.init(cParams);
            obj.createMatrices(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function sp = computeSP(obj,f,g)
            spM = obj.computeSP_M(f,g);
            spK = obj.computeSP_K(f,g);
            sp  = obj.epsilon^2*spK + spM;
        end
        
        function sp = computeSP_M(obj,f,g)
            sp = obj.computeProduct(obj.Msmooth,f,g);
        end
        
        function sp = computeSP_K(obj,f,g)
            sp = obj.computeProduct(obj.Ksmooth,f,g);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.epsilon = cParams.epsilon;
            obj.nVariables = cParams.nVariables;
        end
        
        function createMatrices(obj,cParams)
            obj.mesh = cParams.mesh;
            M = obj.computeMassMatrix();
            K = obj.computeStiffnessMatrix();
            obj.Ksmooth = K;
            obj.Msmooth = M;
        end
        
        function n = computeProduct(obj,K,f,g)
            nx = length(f)/obj.nVariables;
            n = 0;
            for ivar = 1:obj.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                fs = f(i0:iF);
                gs = g(i0:iF);
                n = n + fs'*K*gs;
            end
        end
        
        function M = computeMassMatrix(obj)
            s.type  = 'MassMatrixFun';
            s.mesh  = obj.mesh;
            s.fun   = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATICMASS';
            LHS = LHSintegrator.create(s);
            M = LHS.compute();
        end
    
        function K = computeStiffnessMatrix(obj)
            s.mesh    = obj.mesh;
            s.fValues = zeros(obj.mesh.nnodes, 1);
            f = P1Function(s);
            s.type  = 'StiffnessMatrixFun';
            s.mesh  = obj.mesh;
            s.fun = f;
            LHS = LHSintegrator.create(s);
            K = LHS.compute();
        end
        
    end
    
end
