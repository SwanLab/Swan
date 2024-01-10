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
        materialInterpolator
    end

    methods (Access = public)

        function obj = ComplianceFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function computeFunctionAndGradient(obj,x)
            C = obj.computeMaterial(x);
            u = obj.computeStateVariable(C);            
            J = obj.computeFunction(C,u);
            dC = obj.computeMaterialDerivative(x);   
            dJ = obj.computeGradient(dC,u);
            obj.value    = J; 
            obj.gradient = dJ;
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                 = cParams.mesh;
            obj.filter               = cParams.filter;
            obj.stateProblem         = cParams.stateProblem;
            obj.materialInterpolator = cParams.materialInterpolator;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function J = computeFunction(obj,C,u)
            strain     = obj.computeStateStrain(u);
            stress     = obj.computeStress(C,strain);
            s.mesh     = obj.mesh;
            s.quadType = obj.quadrature.order;
            int        = IntegratorScalarProduct(s);   %create         
            J          = int.compute(strain,stress);
        end

        function u = computeStateVariable(obj,C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u  = obj.stateProblem.uFun;            
        end

        function eu = computeStateStrain(obj,u)
            eu = u.computeSymmetricGradient(obj.quadrature);
            eu = eu.obtainVoigtFormat();
        end

        function stress = computeStress(obj,C,strain)
            Cij = C.evaluate(obj.quadrature.posgp);
            strainV(:,1,:,:) = strain.fValues;
            stress = pagemtimes(Cij,strainV);
            stress = permute(stress, [1 3 4 2]);
            stress = obj.createGaussFunction(stress);
        end

        function g = computeGradient(obj,x)
            dj = obj.computeDJ(x);
            g  = obj.filter.compute(dj,'LINEAR'); 
        end        

        function dj = computeDJ(obj,dC,u)
            dCij         = dC.evaluate(obj.quadrature.posgp);
            eu           = obj.computeStateStrain(u);
            euj(:,1,:,:) = eu.fValues;
            eui(1,:,:,:) = eu.fValues;            
            dStress      = pagemtimes(dCij,euj);
            dj           = pagemtimes(eui,dStress);
            dj           = squeezeParticular(-dj,1);                        
            dj           = obj.createGaussFunction(dj);
        end

        function C = computeMaterial(obj,x)
            mI = obj.materialInterpolator;
            C  = mI.compute(x);
        end

        function dC = computeMaterialDerivative(obj,x)
            mI = obj.materialInterpolator;
            dC = mI.computeDerivative(x);
        end

        function fd = createDiscontinousGaussFunction(obj,f)
            m  = obj.mesh;
            q  = obj.quadrature;
            fd = FGaussDiscontinuousFunction.create(f,m,q);
        end

    end
end