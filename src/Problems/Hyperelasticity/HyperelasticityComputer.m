classdef HyperelasticityComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material
        functional
    end

    properties (Access = private)
        monitor
        updater
    end

    methods (Access = public)

        function obj = HyperelasticityComputer(cParams)
            obj.init(cParams);
            obj.setMonitoring(cParams);
            obj.setOptimizer(cParams)
        end

        function outputData = compute(obj)
            uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            cost = 0;

            nSteps = length(obj.boundaryConditions.bcValues);
            for iStep = 1:nSteps
                [uFun,bc] = obj.preprocess(iStep,nSteps,uFun);
                [uFun,rFun,cost,iterMax] = obj.updater.update(uFun,bc,cost);
                obj.postprocess(iStep,cost,uFun,rFun,iterMax);
            end
            outputData = obj.monitor.data;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material           = cParams.material;
            obj.functional         = cParams.functional;
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay   = cParams.monitoring.set;
            s.shallPrintInfo = cParams.monitoring.printInfo;
            s.shallPrintFile = cParams.monitoring.printFile;
            s.fileNameOut    = cParams.monitoring.fileNameOut;
            obj.monitor = HyperelasticityMonitoring(s);
        end

        function setOptimizer(obj,cParams)
            s.functional = obj.functional;
            s.monitor    = obj.monitor;
            s.tolerance  = cParams.tolerance;
            s.maxIter    = cParams.maxIter;
            obj.updater = DisplacementUpdater(s);
        end

        function [uFun,bc] = preprocess(obj,iStep,nSteps,uFun)
            obj.monitor.printStep(iStep,nSteps)
            bc = obj.boundaryConditions.nextStep();
            uFun.setFValues(obj.updateInitialDisplacement(uFun,bc));
        end

        function u = updateInitialDisplacement(~,uOld,bc)
            restrictedDofs = bc.dirichlet_dofs;
            if isempty(restrictedDofs)
                u = uOld;
            else
                dirich = bc.dirichletFun;
                uVec = reshape(uOld.fValues',[uOld.nDofs 1]);
                dirichVec = reshape(dirich.fValues',[dirich.nDofs 1]);

                uVec(restrictedDofs) = dirichVec(restrictedDofs);
                u = reshape(uVec,[flip(size(uOld.fValues))])';
            end
        end

        function uElas = computeElasticProblem(obj,bc)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material.tensor;
            s.dim = '2D';
            s.boundaryConditions = bc;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            fem.solve();
            uElas = fem.uFun;
        end

        function [rBC,uBC] = postprocess(obj,iStep,cost,uFun,rFun,iterMax)
            [uBC,rBC] = obj.computeNorm(uFun,rFun);
            lambdas   = obj.computeStretches(uFun);
            sigma     = obj.computeCauchyStress(uFun);
            obj.printAndSave(iStep,cost(end),uBC,rBC,uFun,rFun,iterMax);
            obj.monitor.printOutput(iStep,uFun,rFun);
        end

        function [uBC,rBC] = computeNorm(~,uFun,rFun)
            rBC = Norm(rFun,'L2');
            uBC = Norm(uFun,'L2');
        end

        function lambdas = computeStretches(~,u)
            F = Identity(u) + Grad(u);
            C = F'*F;
            lambdas = (Eigen(C).^0.5);
        end

        function sigma = computeCauchyStress(obj,u)
            mu     = obj.material.prop.mu;
            lambda = obj.material.prop.lambda;
            I = Identity(u);
            F = I + Grad(u);
            b = F*F';
            jac = Det(F);
            sigma = Expand(mu./jac,2).*(b-I) + ...
                    Expand(lambda.*(log(jac))./jac,2).*I;
        end

        function printAndSave(obj,step,energy,uBC,rBC,uFun,rFun,iterMax)
            obj.monitor.updateAndRefresh(step,{[energy;uBC],[],[rBC;uBC],[iterMax]});
            obj.saveData(step,energy,uBC,rBC,uFun,rFun,iterMax);
        end

        function saveData(obj,step,energy,uBC,rBC,uFun,rFun,iterMax)
            s.r    = rBC;
            s.rFun = rFun;
            s.u    = uBC;
            s.uFun = uFun;
            s.energy = energy;
            s.numIter = iterMax;
            obj.monitor.saveData(step,s)
        end

    end

end
