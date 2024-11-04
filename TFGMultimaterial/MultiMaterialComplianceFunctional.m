classdef MultiMaterialComplianceFunctional < handle

    properties (Access = private)
        mesh
        filter
        stateProblem
        nMat
        materialInterpolator
        material
    end

    properties (Access = private)
        value0
    end

    methods (Access = public)
        
        function obj = MultiMaterialComplianceFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            obj.material.setDesignVariable(xR);
            C = obj.material.obtainTensor();
            u  = obj.computeStateVariable(C);
            J  = obj.computeFunction(C,u);
            dJc = obj.computeGradient(u,x);
            dJv = [];
            for i = 1:length(dJc)
                dJv = [dJv;dJc{i}.fValues];
            end
            dJ.fValues = dJv;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                 = cParams.mesh;
            obj.filter               = cParams.filter;
            obj.stateProblem         = cParams.stateProblem;
            obj.nMat                 = cParams.nMat;
            obj.materialInterpolator = cParams.materialInterpolator;
            obj.material             = cParams.material;
        end

        function xR = filterFields(obj,x)
            xR = cell(size(x));
            for i = 1:length(x)
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function u = computeStateVariable(obj,C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end

        function J = computeFunction(obj,C,u)
            dCompliance = ElasticEnergyDensity(C,u);
            J           = Integrator.compute(dCompliance,obj.mesh,2);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J = J/obj.value0;
        end

        function dJ = computeGradient(obj,u,x)
            s.mesh           = obj.mesh;
            s.designVariable = x;
            multGrad         = MultimaterialGradientComputer(s);
            TD               = obj.computeTopologicalDerivatives(u);
            dt               = multGrad.compute(TD);
            dJ               = obj.filterFields(dt);
        end

        function dj = computeLocalGradient(obj,u,i,j)
            dC      = obj.material.obtainTensorDerivative(i,j);
            strain  = SymGrad(u);
            dStress = DDP(dC,strain);
            dj      = -0.5.*DDP(strain, dStress);
        end

        function DJ = computeTopologicalDerivatives(obj,u)
            DJ = cell(obj.nMat,obj.nMat);
            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    dj      = obj.computeLocalGradient(u,i,j);
                    DJ{i,j} = dj./obj.value0;
                end
            end
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end