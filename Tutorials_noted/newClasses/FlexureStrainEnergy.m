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
    end

    methods (Access = public)
        function obj = FlexureStrainEnergy(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J, dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filter.compute(xD,2);
            [C,dC] = obj.computeTensorFunctionAndGradient(xR);
            uS = obj.computeStateVariable(C);

            J = obj.computeFunctionValue(C,uS);
            dJ = obj.computeGradient(dC{1},uS);
            dJ = {obj.filter.compute(dJ,2)};

            if ~isempty(obj.value0)
                dJVal = dJ{1}.fValues/obj.value0;
                dJ{1}.setFValues(dJVal);
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.filter       = cParams.filter;
            obj.material     = cParams.material;
            obj.stateProblem = cParams.stateProblem; 

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
        function dJ = computeGradient(dC,uS,uA)
            stateStrain   = SymGrad(uS);
            dStress       = DDP(dC,stateStrain);
            dJ            = DDP(stateStrain,dStress);

            % dE = 1/2 * u * dK * u
            dJ.setFValues(0.5*dJ.fValues)
        end
    end



end


