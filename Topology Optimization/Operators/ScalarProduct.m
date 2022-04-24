classdef ScalarProduct < handle
    
    properties (Access = public)
        epsilon
        Ksmooth
        Msmooth
    end
    
    properties (Access = private)
       nVariables
       mesh
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
            dim = obj.computeDimensions();
            M = obj.computeMassMatrix(dim);
            K = obj.computeStiffnessMatrix(dim);
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

        function dim = computeDimensions(obj)
            s.name = 'x';
            s.mesh = obj.mesh;
            dim = DimensionScalar(s);
        end
        
        function M = computeMassMatrix(obj, dim)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = dim;
            LHS = LHSintegrator.create(s);
            M = LHS.compute();
        end
    
        function K = computeStiffnessMatrix(obj, dim)
            s.type = 'StiffnessMatrix';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = dim;
            LHS = LHSintegrator.create(s);
            K = LHS.compute();
        end
        
    end
    
end
