classdef MaterialInterpolator < handle
    
   properties (Access = protected)
        muFunc
        dmuFunc
        kappaFunc
        dkappaFunc
   end

   properties (Access = protected)
        ndim
        pdim
        matA
        matB
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f   = MaterialInterpolatorFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function m = compute(obj,rho)
            mu      = obj.muFunc(rho.fValues);
            kappa   = obj.kappaFunc(rho.fValues);
            m = obj.createMaterial(mu,kappa);
        end

        function m = computeDerivative(obj,rho)
            dmu       = obj.dmuFunc(rho.fValues);
            dkappa    = obj.dkappaFunc(rho.fValues);         
            m = obj.createMaterial(dmu,dkappa);            
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.pdim = cParams.dim;
            obj.matA = cParams.matA;
            obj.matB = cParams.matB;
            obj.computeNDim(cParams);
        end

        function computeNDim(obj,cParams)
            switch cParams.dim
                case '2D'
                  obj.ndim = 2;
                case '3D'
                  obj.ndim = 3;
            end
        end  
        
        function computeSymbolicInterpolationFunctions(obj)
            [muS,dmuS,kS,dkS] = obj.computeSymbolicMuKappa();
            obj.muFunc        = matlabFunction(muS);
            obj.dmuFunc       = matlabFunction(dmuS);
            obj.kappaFunc     = matlabFunction(kS);
            obj.dkappaFunc    = matlabFunction(dkS);
        end
               
    end

    methods (Access = private)

        function [muS,dmuS,kS,dkS] = computeSymbolicMuKappa(obj)
            [muS,dmuS] = obj.computeMuSymbolicFunctionAndDerivative();
            [kS,dkS]   = obj.computeKappaSymbolicFunctionAndDerivative();
        end

        function m = createMaterial(obj,mu,kappa)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.ndim;
            s.shear   = mu;            
            s.bulk    = kappa;
            m = Material.create(s);   
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeMuSymbolicFunctionAndDerivative(obj)
        computeKappaSymbolicFunctionAndDerivative(obj)
    end
    
end