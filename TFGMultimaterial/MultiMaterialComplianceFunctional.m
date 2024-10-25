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
            xR = obj.filterDesignVariable(xD);
            obj.material.setDesignVariable(xR);
            C = obj.material.obtainTensor();
            u  = obj.computeStateVariable(C);
            J  = obj.computeFunction(C,u);
            dJ = obj.computeGradient(u,x);
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

        function xR = filterDesignVariable(obj,x)
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
            dJval            = pdeprtni(obj.mesh.coord',obj.mesh.connec',dt);
            dJval            = reshape(dJval,[],1);
            dJ.fValues       = dJval;
        end

        function dj = computeLocalGradient(obj,u)
            dC      = obj.material.obtainTensorDerivative();
            strain  = SymGrad(u);
            dStress = DDP(dC,strain);
            dj      = -0.5.*DDP(strain, dStress);
        end

        function DJ = computeTopologicalDerivatives(obj,u)
            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    obj.materialInterpolator.computeFromTo(i,j);
                    dj      = obj.computeLocalGradient(u);
                    djR     = obj.filter.compute(dj,2);
                    derTop  = djR.fValues';
                    TD      = derTop;
                    DJ{i,j} = TD/obj.value0;
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