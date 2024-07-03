classdef ComplianceFunctional < handle

    properties (Access = public)
        bulkValue
        shearValue
    end

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        compliance
        material
    end

    methods (Access = public)
        function obj = ComplianceFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            obj.material.setDesignVariable(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient();
        end

    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.material   = cParams.material;
            obj.compliance = cParams.complainceFromConstitutive;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,'LINEAR');
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj)
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            [u,J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
            dJ     = obj.filter.compute(dJ,'LINEAR');
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J          = obj.computeNonDimensionalValue(J);
            dJ.fValues = obj.computeNonDimensionalValue(dJ.fValues);

            [mu,kappa] = obj.material.obtainShearBulk();
            % Bulk compliance:
            divu  = Divergence(u);
            dbC   = DDP(kappa,DDP(divu,divu));
            bC    = Integrator.compute(dbC,obj.mesh,'QUADRATIC');
            bC    = obj.computeNonDimensionalValue(bC);
            obj.bulkValue = bC;

            % Shear compliance
            e   = AntiVoigt(SymGrad(u));
            D   = Voigt(Deviatoric(e));
            A   = VoigtDeviatorNormMaterial(obj.mesh);
            dsC = DDP(mu,DDP(D,DDP(A,D)));
            sC  = Integrator.compute(dsC,obj.mesh,'QUADRATIC');
            sC  = obj.computeNonDimensionalValue(sC);
            obj.shearValue = sC;
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end