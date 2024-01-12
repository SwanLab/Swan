classdef ComplianceFunctional < handle

    properties (Access = private)
        quadrature
        value0
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

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xR = obj.filterDomainVariable(x);
            C  = obj.computeMaterial(xR);
            u  = obj.computeStateVariable(C);            
            J  = obj.computeFunction(C,u);
            dC = obj.computeMaterialDerivative(xR);   
            dJ = obj.computeGradient(dC,u);
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

        function xR = filterDomainVariable(obj,x)
            xR = obj.filter.compute(x,'LINEAR');
        end

        function C = computeMaterial(obj,x)
            mI = obj.materialInterpolator;
            C  = mI.computeConsitutiveTensor(x);
        end

        function u = computeStateVariable(obj,C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end

        function J = computeFunction(obj,C,u)
            strain = obj.computeStateStrain(u);
            stress = obj.computeStress(C,strain);
            int    = Integrator.create('ScalarProduct',obj.mesh,obj.quadrature.order);
            J      = int.compute(strain,stress);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = obj.computeNonDimensionalValue(J);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
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

        function dC = computeMaterialDerivative(obj,x)
            mI = obj.materialInterpolator;
            dC = mI.computeConsitutiveTensorDerivative(x);
        end

        function g = computeGradient(obj,dC,u)
            dj = obj.computeDJ(dC,u);
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
            dj           = obj.computeNonDimensionalValue(dj);
            dj           = obj.createGaussFunction(dj);
        end

        function fd = createGaussFunction(obj,f)
            m  = obj.mesh;
            q  = obj.quadrature;
            fd = FGaussDiscontinuousFunction.create(f,m,q);
        end
    end
end