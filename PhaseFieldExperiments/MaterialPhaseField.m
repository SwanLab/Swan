classdef MaterialPhaseField < Material

    properties (Access = private)
        fun
        phi
    end

    properties (Access = private)
        isoMat
        mesh
        matInterpolation
        Gc
    end

    methods (Access = public)

        function obj = setMaterial(obj,phi,type)
            obj.phi = phi;
            switch type
                case 'Isotropic'
                    obj.fun = @(phi) 1;
                case 'Interpolated'
                    obj.fun = obj.matInterpolation.fun;
                case 'Jacobian'
                    obj.fun = obj.matInterpolation.dfun;
                case 'Hessian'
                    obj.fun = obj.matInterpolation.ddfun;
            end
        end

        function obj = MaterialPhaseField(cParams)
            obj.init(cParams)
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

        function init(obj,cParams)
            obj.isoMat = cParams.isoMat;
            obj.mesh = cParams.mesh;
            obj.matInterpolation = cParams.materialInterpolation;
            obj.Gc = cParams.Gc;
        end

    end

    methods (Access = private)

        function [mu,l] = computeMatParams(obj,xV)
            gV = obj.evaluateDegradationFun(xV);
            [mu, l] = obj.computeMuAndLambda(gV,xV);
        end

        function gV = evaluateDegradationFun(obj,xV)
            s.mesh = obj.mesh;
            s.handleFunction = obj.fun;
            s.l2function = obj.phi;
            g = CompositionFunction(s);
            gV = g.evaluate(xV);
            gV = squeezeParticular(gV,1);
        end

        function [mu, l] = computeMuAndLambda(obj,gV,xV)
            [muV,kV] = obj.isoMat.computeShearAndBulk(xV);
            mu = gV.*muV;
            k = gV.*kV;
            l = obj.isoMat.computeLambdaFromShearAndBulk(mu,k,obj.ndim);
        end
    end

end