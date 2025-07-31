classdef HyperelasticityComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
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
                obj.monitor.printStep(iStep,nSteps)
                [uFun,bc] = obj.updateBoundaryConditions(uFun);
                [uFun,rFun,cost,iterMax] = obj.updater.update(uFun,bc,cost);
                [uBC,rBC] = obj.postprocess(uFun,rFun);
                obj.printAndSave(iStep,cost(end),uBC,rBC,uFun,rFun,iterMax);
                obj.monitor.printOutput(iStep,uFun,rFun);
            end
            outputData = obj.monitor.data;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
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
            s.functional  = obj.functional;
            s.monitor     = obj.monitor;
            s.tolerance   = cParams.tolerance;
            s.maxIter     = cParams.maxIter;
            obj.updater = DisplacementUpdater(s);
        end

        function [u,bc] = updateBoundaryConditions(obj,u)
            bc = obj.boundaryConditions.nextStep();
            u.setFValues(obj.updateInitialDisplacement(u,bc));
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

        function [rBC,uBC] = postprocess(~,uFun,rFun)
            rBC = Norm(rFun,'L2');
            uBC = Norm(uFun,'L2');
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

        % 
        % function u = computeInitialElasticGuess(obj, bc)
        %     s.mesh = obj.mesh;
        %     s.scale = 'MACRO';
        %     s.material = obj.material;
        %     s.dim = '2D';
        %     s.boundaryConditions = bc;
        %     s.solverType = 'REDUCED';
        %     s.solverMode = 'DISP';
        %     s.solverCase = 'DIRECT';
        %     fem = ElasticProblem(s);
        %     fem.solve();
        %     %             u = obj.reshapeToVector(fem.uFun.fValues);
        %     u = fem.uFun;
        % end
        % 
        % Other useful stuff later on
        % function lambdas = computeStretches(obj)
        %     GradU2 = ActualGrad(obj.uFun);
        %     Id = Identity(obj.uFun);
        %     F = GradU2 + Id;
        %     C = F'*F;
        %     lambdas = (Eigen(C).^0.5);
        % end
        % 
        % function sigma = computeCauchyStress(obj)
        %     mu = obj.material.mu;
        %     lambda = obj.material.lambda;
        %     GradU2 = ActualGrad(obj.uFun);
        %     Id = Identity(obj.uFun);
        %     F = GradU2 + Id;
        %     b = F*F';
        %     jac = Det(F);
        %     sigma = mu.*(b-Id)./jac + lambda.*(log(jac)).*Id./jac;
        %     sigma.ndimf = [2 2];
        % end


    end

end
