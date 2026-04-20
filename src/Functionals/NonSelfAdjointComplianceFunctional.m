classdef NonSelfAdjointComplianceFunctional < handle

    properties (Access = public)
        cost_in_norm
        cost_out_norm
        value0
    end

        
 
    properties (Access = private)
       % value0
        integralK_value
    end

    properties (Access = private)
        quadrature
        adjointProblem
    end

    properties (Access = private)
        mesh
        filter
        material
        stateProblem
    end

    methods (Access = public)
        function obj = NonSelfAdjointComplianceFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD     = x.obtainDomainFunction();
            xR     = obj.filterDesignVariable(xD);
            [C,dC] = obj.computeTensorFunctionAndGradient(xR);
            uS     = obj.computeStateVariable(C);
            uA    = obj.computeAdjointVariable(C);
            J      = obj.computeFunctionValue(C,uS,uA);
            dJ     = obj.computeGradient(dC{1},uS,uA);
            dJ     = {obj.filter.compute(dJ,2)};
            dJVal  = dJ{1}.fValues/obj.value0;
            dJ{1}.setFValues(dJVal);

        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.filter         = cParams.filter;
            obj.material       = cParams.material;
            obj.stateProblem   = cParams.stateProblem;
            obj.adjointProblem = cParams.adjointProblem;
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh, 2);
            obj.quadrature = quad;
        end

        function xR = filterDesignVariable(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function [C,dC] = computeTensorFunctionAndGradient(obj,xR)
            obj.material.setDesignVariable(xR);
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();
        end

        function u = computeStateVariable(obj,C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end

        function u = computeAdjointVariable(obj,C)
            obj.adjointProblem.updateMaterial(C);
            obj.adjointProblem.solve();
            u = obj.adjointProblem.uFun;
        end

        function J = computeFunctionValue(obj,C,uS,uA)
            stateStrain   = SymGrad(uS);
            adjointStrain = SymGrad(uA);
            stressA       = DDP(C,adjointStrain);
            dCompliance   = DDP(stateStrain,stressA);
            J             = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);

            t = obj.adjointProblem.boundaryConditions.tractionFun; % loads of the adjoint problem
            rhs_in  = t(1).computeRHS(obj.adjointProblem.uFun); % input loads 
            rhs_out = t(2).computeRHS(obj.adjointProblem.uFun); % output loads
            t_in_reshaped  = reshape(rhs_in, 2, [])';
            t_out_reshaped = reshape(rhs_out, 2, [])';

            t_in_lagrangian = LagrangianFunction.create(obj.mesh,2,'P1');
            t_out_lagrangian = LagrangianFunction.create(obj.mesh,2,'P1');

            t_in_lagrangian.setFValues(t_in_reshaped);
            t_out_lagrangian.setFValues(t_out_reshaped);

            % u_vec = reshape(uS.fValues',[],1); % to undo the reshape done in computeDisplacements

            cost_in = sum(Integrator.compute(DP(t_in_lagrangian,uS),obj.mesh,2));
            cost_out = sum(Integrator.compute(DP(t_out_lagrangian,uS),obj.mesh,2));

            int_k_in = sum(Integrator.compute(abs(t_in_lagrangian),obj.mesh,2));
            int_k_out = sum(Integrator.compute(abs(t_out_lagrangian),obj.mesh,2));  

            % int_k_in = Integrator.compute(abs(rhs_in),obj.mesh,2);
            % int_k_out = Integrator.compute(abs(rhs_out),obj.mesh,2);  


            %fprintf('Iter: worker computed In: %f, Out: %f\n', obj.cost_in_norm, obj.cost_out_norm);

            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = J/obj.value0;

            value0vec = [obj.value0, 0];
            obj.cost_in_norm = cost_in / value0vec;
            obj.cost_out_norm = cost_out / value0vec;
        end
    end

    methods (Static, Access = private)
        function dJ = computeGradient(dC,uS,uA)
            stateStrain   = SymGrad(uS);
            adjointStrain = SymGrad(uA);
            dStressA      = DDP(dC,adjointStrain);
            dJ            = -DDP(stateStrain,dStressA);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'NonSelfAdjCompliance';
        end
    end
end