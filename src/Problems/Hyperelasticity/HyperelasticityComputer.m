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
        maxIter
        stopWhenNotConverged
        usePredictor
        checkOverlaps
        overlapCounter
        interpolationOrder
        maxBisections
    end

    methods (Access = public)

        function obj = HyperelasticityComputer(cParams)
            obj.init(cParams);
            obj.setMonitoring(cParams);
            obj.setOptimizer(cParams)
        end

        function outputData = compute(obj)
            uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, obj.interpolationOrder);
            cost = 0;

            nSteps = length(obj.boundaryConditions.bcValues);
            for iStep = 1:nSteps
                obj.monitor.printStep(iStep,nSteps)
                [uFun,rFun,cost,iterMax] = obj.solveStep(iStep,uFun,cost);
                if obj.stopWhenNotConverged && obj.hasFailed(uFun,iterMax)
                    break
                end
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
            obj.maxIter            = cParams.maxIter;
            obj.stopWhenNotConverged = obj.readOption(cParams,'stopWhenNotConverged',false);
            obj.usePredictor         = obj.readOption(cParams,'usePredictor',false);
            obj.checkOverlaps        = obj.readOption(cParams,'checkOverlaps',false);
            obj.interpolationOrder   = obj.readOption(cParams,'interpolationOrder','P1');
            obj.maxBisections        = obj.readOption(cParams,'maxBisections',0);
            obj.createOverlapCounter();
        end

        function createOverlapCounter(obj)
            if obj.checkOverlaps
                s.mesh    = obj.mesh;
                s.penalty = 0;
                s.depth   = 0;
                obj.overlapCounter = SelfContactFunctional(s);
            end
        end

        function [uFun,rFun,cost,iter] = solveStep(obj,iStep,uFun,cost)
            if obj.maxBisections == 0
                bc = obj.boundaryConditions.nextStep();
                [uFun,rFun,cost,iter] = obj.solveFor(bc,uFun,cost);
                return
            end
            values = obj.boundaryConditions.bcValues;
            target = values(iStep);
            start  = values(max(iStep-1,1));
            [uFun,rFun,cost,iter] = obj.solveBetween(start,target,uFun,cost,0);
        end

        function [uFun,rFun,cost,iter] = solveBetween(obj,start,target,uFun,cost,depth)
            uStart = copy(uFun);
            bc = obj.boundaryConditions.computeConditions(target);
            [uFun,rFun,cost,iter] = obj.solveFor(bc,uFun,cost);
            if ~obj.hasFailed(uFun,iter) || depth == obj.maxBisections
                return
            end
            middle = (start + target)/2;
            [uFun,rFun,cost,iter] = obj.solveBetween(start,middle,uStart,cost,depth+1);
            if ~obj.hasFailed(uFun,iter)
                [uFun,rFun,cost,iter] = obj.solveBetween(middle,target,uFun,cost,depth+1);
            end
        end

        function [uFun,rFun,cost,iter] = solveFor(obj,bc,uFun,cost)
            uFun = obj.preprocess(bc,uFun);
            [uFun,rFun,cost,iter] = obj.updater.update(uFun,bc,cost);
            [uFun,rFun,cost,iter] = obj.enforceParameters(uFun,rFun,cost,iter,bc);
        end

        function [uFun,rFun,cost,iter] = enforceParameters(obj,uFun,rFun,cost,iter,bc)
            if ~obj.hasParameters()
                return
            end
            for k = 1:10
                if iter >= obj.maxIter || obj.functional.updateParameters(uFun)
                    return
                end
                [uFun,rFun,cost,iter] = obj.updater.update(uFun,bc,cost);
            end
        end

        function itHas = hasParameters(obj)
            itHas = ismethod(obj.functional,'updateParameters');
        end

        function itHas = hasFailed(obj,uFun,iter)
            u = uFun.fValues;
            itHas = iter >= obj.maxIter || ~isreal(u) || any(~isfinite(u(:)));
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
            s.lineSearch = isfield(cParams,'lineSearch') && cParams.lineSearch;
            if s.lineSearch
                s.residualTolerance = cParams.residualTolerance;
                s.viscosity         = cParams.viscosity;
            end
            obj.updater = DisplacementUpdater(s);
        end

        function uFun = preprocess(obj,bc,uFun)
            if obj.hasParameters()
                obj.functional.resetParameters();
            end
            uVec = obj.toDofVector(uFun);
            du   = obj.computeDirichletIncrement(uVec,bc);
            if obj.usePredictor
                du = obj.addPredictedFreeIncrement(uFun,uVec,du,bc);
            end
            uFun.setFValues(obj.toNodalValues(uVec+du,uFun));
        end

        function du = addPredictedFreeIncrement(obj,uFun,uVec,du,bc)
            K     = obj.functional.computeHessian(uFun);
            free  = bc.free_dofs;
            dDofs = bc.dirichlet_dofs;
            duFree = DisplacementUpdater.solveWithPositiveShift(K(free,free),-K(free,dDofs)*du(dDofs));
            eBest  = obj.computeEnergyOf(uVec+du,uFun,bc);
            for alpha = [1 0.5 0.25 0.125]
                candidate       = du;
                candidate(free) = alpha*duFree;
                e = obj.computeEnergyOf(uVec+candidate,uFun,bc);
                if obj.isValidEnergy(e) && (~obj.isValidEnergy(eBest) || e < eBest)
                    du    = candidate;
                    eBest = e;
                end
            end
        end

        function e = computeEnergyOf(obj,uVec,uFun,bc)
            trial = copy(uFun);
            trial.setFValues(obj.toNodalValues(uVec,uFun));
            e = obj.functional.computeCost(trial,bc);
        end

        function du = computeDirichletIncrement(obj,uVec,bc)
            dDofs     = bc.dirichlet_dofs;
            dVec      = obj.toDofVector(bc.dirichletFun);
            du        = zeros(size(uVec));
            du(dDofs) = dVec(dDofs) - uVec(dDofs);
        end

        function itIs = isValidEnergy(~,e)
            itIs = isreal(e) && isfinite(e);
        end

        function v = toDofVector(~,fun)
            v = reshape(fun.fValues',[fun.nDofs 1]);
        end

        function f = toNodalValues(~,v,fun)
            f = reshape(v,[flip(size(fun.fValues))])';
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
            [uBC,rBC] = obj.computeNorm(iStep,uFun,rFun);
            lambdas   = obj.computeStretches(uFun);
            sigma     = obj.computeCauchyStress(uFun);
            obj.printAndSave(iStep,cost(end),uBC,rBC,uFun,rFun,iterMax);
            obj.monitor.printOutput(iStep,uFun,rFun);
        end

        function [uBC,rBC] = computeNorm(obj,step,uFun,rFun)
            UpSide = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            if ismember(obj.boundaryConditions.type, ["ForceTractionY", "ForceTractionYClamped"])
                uBC = norm(mean(uFun.fValues(nodes(isInUp),2)));
                rBC = obj.boundaryConditions.bcValues(step);
            elseif ismember(obj.boundaryConditions.type, ["DisplacementTractionY","DisplacementTractionYClamped"])
                rBC = norm(sum(rFun.fValues(nodes(isInUp),2)));
                uBC = obj.boundaryConditions.bcValues(step);
            else
                rBC = Norm(rFun,'L2');
                uBC = Norm(uFun,'L2');
            end
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
            s.uFun = copy(uFun);
            s.energy = energy;
            if obj.checkOverlaps
                s.overlaps = obj.overlapCounter.countOverlaps(uFun);
            end
            s.numIter = iterMax;
            obj.monitor.saveData(step,s)
        end

    end

    methods (Static, Access = private)

        function value = readOption(cParams,name,default)
            value = default;
            if isfield(cParams,name)
                value = cParams.(name);
            end
        end

    end

end
