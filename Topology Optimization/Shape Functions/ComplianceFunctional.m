classdef ComplianceFunctional < handle

    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        quadrature
    end

    properties (Access = private)
        mesh
        filter
        stateProblem
        adjointProblem
        materialInterpolator
    end

    methods (Access = public)

        function obj = ComplianceFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function computeFunctionAndGradient(obj,x)
            obj.computeFunction(x);
            obj.computeGradient(x);
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                 = cParams.mesh;
            obj.filter               = cParams.filter;
            obj.stateProblem         = cParams.stateProblem;
            obj.adjointProblem       = obj.stateProblem; 
            obj.materialInterpolator = cParams.materialInterpolator;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function computeFunction(obj,x)
            s.mesh     = obj.mesh;
            s.quadType = obj.quadrature.order;
            int        = IntegratorScalarProduct(s);
            C          = obj.computeMaterial(x);
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            strain     = obj.computeStateStrain();
            stress     = obj.computeStress(C,strain);
            J          = int.compute(strain,stress);
            obj.value  = J;
        end

        function strain = computeStateStrain(obj)
            u      = obj.stateProblem.uFun;
            strain = u.computeSymmetricGradient(obj.quadrature);
            strain = strain.obtainVoigtFormat();
        end

        function strain = computeAdjointStrain(obj)
            p      = obj.adjointProblem.uFun;
            strain = p.computeSymmetricGradient(obj.quadrature);
            strain = strain.obtainVoigtFormat();
        end

        function stress = computeStress(obj,C,strain)
            Cij              = C.evaluate(obj.quadrature.posgp);
            strainV(:,1,:,:) = strain.fValues;
            stressV          = pagemtimes(Cij,strainV);
            stressV          = permute(stressV, [1 3 4 2]);
            m                = obj.mesh;
            q                = obj.quadrature;
            stress           = FGaussDiscontinuousFunction.create(stressV,m,q);
        end

        function computeGradient(obj,x)
            dj = obj.computeDJ(x);
            g  = obj.filter.compute(dj,'LINEAR'); 
            obj.gradient = g;
        end

        function dj = computeDJ(obj,x)
            dC           = obj.computeMaterialDerivative(x);   
            dCij         = dC.evaluate(obj.quadrature.posgp);
            eu           = obj.computeStateStrain();
            ep           = obj.computeAdjointStrain();
            epi(1,:,:,:) = ep.fValues;
            euj(:,1,:,:) = eu.fValues;
            dStress      = pagemtimes(dCij,euj);
            dj           = pagemtimes(epi,dStress);
            dj           = squeezeParticular(-dj,1);                        
            m            = obj.mesh;
            q            = obj.quadrature;
            dj           = FGaussDiscontinuousFunction.create(dj,m,q);
        end

        function C = computeMaterial(obj,x)
            mI = obj.materialInterpolator;
            C  = mI.compute(x);
        end

        function dC = computeMaterialDerivative(obj,x)
            mI = obj.materialInterpolator;
            dC = mI.computeDerivative(x);
        end
    end
end