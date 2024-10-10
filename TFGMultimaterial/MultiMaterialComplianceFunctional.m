classdef MultiMaterialComplianceFunctional < handle

    properties (Access = private)
        mesh
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
            obj.material.setDesignVariable(x);
            C = obj.material.obtainTensor();
            u  = obj.computeStateVariable(C);
            J  = obj.computeFunction(C,u);
            dJ = obj.computeGradient(u,x);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                 = cParams.mesh;
            obj.stateProblem         = cParams.stateProblem;
            obj.nMat                 = cParams.nMat;
            obj.materialInterpolator = cParams.materialInterpolator;
            obj.material             = cParams.material;
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
            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    obj.materialInterpolator.computeFromTo(i,j);
                    dj      = obj.computeLocalGradient(u);
                    derTop  = squeezeParticular(dj.evaluate([0;0]),2);
                    TD      = derTop;
                    DJ{i,j} = TD/obj.value0;
                end
            end
            dJ = obj.smoothGradient(DJ,x);
            dJ = reshape(dJ,[],1);
        end

        function dj = computeLocalGradient(obj,u)
            dC      = obj.material.obtainTensorDerivative();
            strain  = SymGrad(u);
            dStress = DDP(dC,strain);
            dj      = -0.5.*DDP(strain, dStress);
        end
        
        function dJ = smoothGradient(obj,TD,x)
            s.mesh           = obj.mesh;
            s.designVariable = x;
            multGrad         = MultimaterialGradientComputer(s);
            dt               = multGrad.compute(TD);
            t                = obj.mesh.connec';
            p                = obj.mesh.coord';
            dJ               = pdeprtni(p,t,dt);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end