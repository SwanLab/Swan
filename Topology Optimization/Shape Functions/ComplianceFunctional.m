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
            xR = obj.filterFields(xD);
            obj.material.setDesignVariable(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient(x);
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

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj,x)
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            dC  = ChainRule.compute(x,dC);
            [u,J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
            dJ     = obj.filterFields(dJ);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J  = obj.computeNonDimensionalValue(J);
            dJ = obj.computeNonDimensionalGradient(dJ);

            [mu,kappa] = obj.material.obtainShearBulk();
            % Bulk compliance:
            divu  = Divergence(u);
            dbC   = DDP(kappa,DDP(divu,divu));
            bC    = Integrator.compute(dbC,obj.mesh,2);
            bC    = obj.computeNonDimensionalValue(bC);
            obj.bulkValue = bC;

            % Shear compliance
            e   = AntiVoigt(SymGrad(u));
            D   = Voigt(Deviatoric(e));
            A   = VoigtDeviatorNormMaterial(obj.mesh);
            dsC = DDP(mu,DDP(D,DDP(A,D)));
            sC  = Integrator.compute(dsC,obj.mesh,2);
            sC  = obj.computeNonDimensionalValue(sC);
            obj.shearValue = sC;
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end

        function dx = computeNonDimensionalGradient(obj,dx)
            refX = obj.value0;
            for i = 1:length(dx)
                dx{i}.setFValues(dx{i}.fValues/refX);
            end
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end