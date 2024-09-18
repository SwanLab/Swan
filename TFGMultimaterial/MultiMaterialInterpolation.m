classdef MultiMaterialInterpolation < handle
    
    % PENDENT DE FER SERVIR

   properties (Access = private)
        muFunc
        dmuFunc
        kappaFunc
        dkappaFunc
   end

   properties (Access = private)
        mat
        C0
   end

    methods (Access = public)
        function obj = MultiMaterialInterpolation(cParams)
            obj.init(cParams)
        end

        function [mu,kappa] = computeConsitutiveTensor(obj,chi)
            tgamma     = (obj.mat.E/obj.mat.E(1))*chi;
            Ceff       = obj.C0*tgamma;
            lambdaVals = Ceff(6,:);
            muVals     = Ceff(4,:);

            s.order     = 'P0';
            s.fValues   = lambdaVals';
            s.mesh      = chi.mesh;
            lambda = LagrangianFunction(s);

            s.order   = 'P0';
            s.fValues = muVals';
            s.mesh    = chi.mesh;
            mu        = LagrangianFunction(s);
            
            N     = chi.mesh.ndim;
            kappa = lambda + 2.*mu/N;
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,rho)
            dmu      = CompositionFunction.create(obj.dmuFunc,rho);
            dkappa   = CompositionFunction.create(obj.dkappaFunc,rho);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mat = cParams.mat;
            obj.C0  = cParams.C0;
        end
    end
end