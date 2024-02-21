classdef MaterialPhaseField < IsotropicElasticMaterial

    properties (Access = private)
        fun
        phi
    end

    properties (Access = private)
        mesh
        matInterpolation
        Gc
    end

    methods (Access = public)

        function obj = MaterialPhaseField(cParams)
            obj = obj@IsotropicElasticMaterial(cParams);
            obj.initPhaseField(cParams)
        end

        function C = evaluateIsotropicMaterial(obj,xV)
            [mu,kappa] = obj.computeShearAndBulk(xV);
            lambda = obj.computeLambdaFromShearAndBulk(mu,kappa,obj.ndim);
            C = obj.evaluateMaterial(mu,lambda);
        end

        function C = evaluate(obj,phi,xV,type)
            switch type
                case 
        function C = evaluateInterpolatedMaterial(obj,phi,xV)
            fun = obj.matInterpolation.fun;
        end

        function C = evaluateFirstDerivativeInterpolatedMaterial(obj,phi,xV)
            dfun = obj.matInterpolation.dfun;
        end

        function C = evaluateSecondDerivativeInterpolatedMaterial(obj,phi,xV)
            ddfun = obj.matInterpolation.ddfun;
        end

    end

    methods (Access = private)

        function initPhaseField(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.matInterpolation = cParams.materialInterpolation;
            obj.Gc = cParams.Gc;
        end

        function C = evaluate(obj,xV)
            [mu,l] = obj.computeMatParams(xV);
            nPoints = size(mu,1);                        
            nElem = size(mu,2);
            nStre = 3;
            C = zeros(nStre,nStre,nPoints,nElem);
            C(1,1,:,:)= 2*mu+l;
            C(1,2,:,:)= l;
            C(2,1,:,:)= l;
            C(2,2,:,:)= 2*mu+l;
            C(3,3,:,:)= mu;
        end

    end

    methods (Access = private)

        function [mu,l] = computeMatParams(obj,xV)
            gV = obj.evaluateDegradationFun(obj.fun,obj.phi,xV);
            [mu, l] = obj.computeMuAndLambda(gV,xV);
        end

        function gV = evaluateDegradationFun(obj,xV)
            s.mesh = obj.mesh;
            s.handleFunction = obj.fun;
            s.l2function = obj.phi;
            g = CompositionFunction(s);
            gV = g.evaluate(xV);

            %val = permute(val,[1 3 2]);
            %val = squeezeParticular(val,1);
        end

        function [mu, l] = computeMuAndLambda(obj,gV,xV)
            [muV,kV] = obj.computeShearAndBulk(xV);
            mu = gV.*muV;
            k = gV.*kV;
            l = obj.computeLambdaFromShearAndBulk(mu,k,obj.ndim);
        end
    end

end