classdef MaterialPhaseField < handle

    properties (Access = public)
        material
        Gc
        %fc
    end

    properties (Access = private)
        mesh
        matInterpolation
        E
        nu
    end

    properties (Access = private)
        muVal
        mu
        kappaVal
        kappa
    end


    methods (Access = public)

        function obj = MaterialPhaseField(cParams)
            obj.init(cParams)
        end

        function computeIsotropicMaterial(obj,quad)
            obj.computeMatIsoParams(quad);
            obj.computeMaterial();
        end

        function computeInterpolatedMaterial(obj,phi,quad)
            obj.computeMatIntParams(phi,quad);
            obj.computeMaterial();
        end

        function computeFirstDerivativeInterpolatedMaterial(obj,phi,quad)
            obj.computeDMatIntParams(phi,quad);
            obj.computeMaterial();
        end

        function computeSecondDerivativeInterpolatedMaterial(obj,phi,quad)
            obj.computeDDMatIntParams(phi,quad);
            obj.computeMaterial();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.matInterpolation = cParams.materialInterpolation;
            obj.E = cParams.E;
            obj.nu = cParams.nu;
            obj.Gc = cParams.Gc;
            %obj.fc = cParams.fc;

            obj.muVal = obj.computeMuFromYoungAndNu();
            obj.kappaVal = obj.computeKappaFromYoungAndNu();
        end

        function computeMatIsoParams(obj,quad)
            obj.mu = obj.muVal*ones(obj.mesh.nelem,quad.ngaus);
            obj.kappa = obj.kappaVal*ones(obj.mesh.nelem,quad.ngaus);
        end

        function computeMatIntParams(obj,phi,quad)
            g = obj.createDegradationFunction(phi,quad);
            obj.kappa = g.*obj.kappaVal;
            obj.mu = g.*obj.muVal;
        end

        function computeDMatIntParams(obj,phi,quad)
            g = obj.createFirstDerivativeDegradationFunction(phi,quad);
            obj.kappa = g.*obj.kappaVal;
            obj.mu = g.*obj.muVal;
        end

        function computeDDMatIntParams(obj,phi,quad)
            g = obj.createSecondDerivativeDegradationFunction(phi,quad);
            obj.kappa = g.*obj.kappaVal;
            obj.mu = g.*obj.muVal;
        end

        function computeMaterial(obj)
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = obj.kappa;
            s.mu    = obj.mu;
            mat = Material.create(s);
            mat.compute(s);

            obj.material = mat;
        end

    end

    methods (Access = private)

        function kappa = computeKappaFromYoungAndNu(obj)
            kappa = obj.E/(2*(1-obj.nu));
        end

        function mu = computeMuFromYoungAndNu(obj)
            mu = obj.E./(2*(1+obj.nu));
        end

        function gVal = createDegradationFunction(obj,phi,quad)
            s.mesh = obj.mesh;
            s.handleFunction = obj.matInterpolation.fun;
            s.l2function = phi;
            g = CompositionFunction(s);
            gVal = obj.computeGaussPointValue(g,quad);
        end

        function dgVal = createFirstDerivativeDegradationFunction(obj,phi,quad)
            s.mesh = obj.mesh;
            s.handleFunction = obj.matInterpolation.dfun;
            s.l2function = phi;
            dg = CompositionFunction(s);
            dgVal = obj.computeGaussPointValue(dg,quad);
        end

        function ddgVal = createSecondDerivativeDegradationFunction(obj,phi,quad)
            s.mesh = obj.mesh;
            s.handleFunction = obj.matInterpolation.ddfun;
            s.l2function = phi;
            ddg = CompositionFunction(s);
            ddgVal = obj.computeGaussPointValue(ddg,quad);
        end

        function val = computeGaussPointValue(obj,fun,quad)
            val = fun.evaluate(quad.posgp);
            val = permute(val,[1 3 2]);
            val = squeezeParticular(val,1);
        end

    end


end