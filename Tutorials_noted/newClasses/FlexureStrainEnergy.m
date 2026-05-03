classdef FlexureStrainEnergy < handle
    % Class to compute the strain energy and its gradient for a given
    % degree

    properties (Access = public)
        value0;
    end

    properties (Access = private)
        mesh
        filter
        material
        stateProblem
        quadrature
        gradientFilter
    end

    methods (Access = public)
        function obj = FlexureStrainEnergy(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J, dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();

            if iscell(xD)
                xD_single = xD{1};
            else
                xD_single = xD;
            end

            x_filtered = obj.filter.compute(xD_single,2);
            xR = {x_filtered};

            [C,dC] = obj.computeTensorFunctionAndGradient(xR);
            uS = obj.computeStateVariable(C);

            J = obj.computeFunctionValue(C,uS);
            dJ = obj.computeGradient(dC{1},uS);
            dJ = obj.gradientFilter.compute(dJ,2);
            dJval = 0.5*dJ.fValues;

            if ~isempty(obj.value0)
                dJval = dJval/obj.value0;
            end

            dJ.setFValues(dJval);
            dJ = {dJ};
        end

        function title = getTitleToPlot(obj)
            title = 'Strain Energy (DOC)';
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.filter       = cParams.filter;
            obj.material     = cParams.material;
            obj.stateProblem = cParams.stateProblem; 
            obj.gradientFilter = cParams.gradientFilter;

            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh, 2);
            obj.quadrature = quad;
        end

        function xR = filterDesignVariable(obj, x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function [C, dC] = computeTensorFunctionAndGradient(obj, xR)
            obj.material.setDesignVariable(xR);
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();
        end
        
        function u = computeStateVariable(obj, C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end
        
        function J = computeFunctionValue(obj, C, uS)
            stateStrain   = SymGrad(uS);
            stress        = DDP(C, stateStrain);
            dCompliance = DDP(stateStrain, stress);
            
            % E = 1/2 * u * K * u
            J = 0.5 * Integrator.compute(dCompliance, obj.mesh, obj.quadrature.order);
            
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = J / obj.value0;
        end

    end

    methods (Static, Access = private)
        function dJ = computeGradient(dC,uS)
            stateStrain   = SymGrad(uS);
            dStress       = DDP(dC,stateStrain);
            dJ            = DDP(stateStrain,dStress);
        end
    end



end


