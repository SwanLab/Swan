classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        internalDamageVariable
        functional
        tolerance
        maxIter
    end

    properties (Access = private)
        monitor
        optimizer
        updater
    end

    methods (Access = public)
        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
            obj.setMonitoring(cParams);
            obj.setOptimizer(cParams);
        end

        function outputData = compute(obj)
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            r = obj.internalDamageVariable;
            cost = 0;

            nSteps = length(obj.boundaryConditions.bcValues);
            for iStep = 1:nSteps
                [u,bc] = obj.preprocess(iStep,nSteps,u);
                [u,r,F,cost,iterMax] = obj.updater.update(u,r,bc,cost);
                obj.postprocess(iStep,u,r,F,cost,iterMax)
            end
            outputData = obj.monitor.data;
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.internalDamageVariable = cParams.internalDamageVariable;
            obj.functional         = cParams.functional;
            obj.tolerance          = cParams.tolerance;
            obj.maxIter            = cParams.maxIter;
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay   = cParams.monitoring.set;
            s.shallPrintInfo = cParams.monitoring.print;
            s.fun            = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.monitor = ContinuumDamageMonitoring(s);
        end

        function setOptimizer(obj,cParams)
            s.functional = obj.functional;
            s.monitor    = obj.monitor;
            s.tolerance  = cParams.tolerance;
            s.maxIter    = cParams.maxIter;
            obj.updater  = ContinuumDamageDisplacementUpdater(s);
        end

        function [u,bc] = preprocess(obj,iStep,nSteps,u)
            obj.monitor.printStep(iStep,nSteps)
            bc = obj.boundaryConditions.nextStep();
            u  = obj.computeInitialDisplacement(u,bc);
        end

        function u = computeInitialDisplacement(obj,u,bc)
            mat = obj.functional.getMaterial();
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = mat.obtainNonDamagedTensor;
            s.dim = '2D';
            s.boundaryConditions = bc;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            fem.solve();
            u.setFValues(fem.uFun.fValues);
        end

        function postprocess(obj,iStep,uFun,r,F,cost,iterMax)
            [fVal,uVal] = obj.computeTotalReaction(iStep,F,uFun);
            [dmgFun0,dmgFun1,qFun,rFun] = obj.computeDamageVariables(uFun,r);
            obj.printAndSave(iStep,uFun,dmgFun0,dmgFun1,qFun,rFun,uVal,fVal,cost(end),iterMax);
        end

        function [fVal,uVal] = computeTotalReaction(obj,step,F,u)
            % Only reactions in Upper Side -> Change for DisplacementX
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            if obj.boundaryConditions.type == "forceTraction"
                uVal = norm(mean(u.fValues(nodes(isInUp),2)));
                fVal = obj.boundaryConditions.bcValues(step);
            else
                ReactX = sum(F(2*nodes(isInUp)-1));
                ReactY = sum(F(2*nodes(isInUp)));
                fVal = sqrt(ReactX^2+ReactY^2);
                uVal = obj.boundaryConditions.bcValues(step);
            end
        end

        function [dmgFunP0,dmgFunP1,qFun,rFun] = computeDamageVariables(obj,u,r)
            dmg = obj.functional.getDamage(r);
            dmgFunP0 = dmg.project('P0');
            dmgFunP1 = dmg.project('P1');

            qFun   = obj.functional.getHardening(r).project('P0');
            rFun   = obj.internalDamageVariable.r.project('P0');
        end

        function printAndSave(obj,step,uFun,dmgFunP0,dmgFunP1,qFun,rFun,uVal,fVal,energy,iterMax)
            dmgMax = max(dmgFunP0.fValues); 
            qMax = max(qFun.fValues);
            rMax = max(rFun.fValues);
            obj.monitor.updateAndRefresh(step,{[fVal;uVal],[dmgMax;uVal],[qMax,rMax],...
                                               [],[dmgFunP1.fValues],[iterMax]});
            obj.saveData(step,uFun,dmgFunP0,qMax,rMax,uVal,fVal,energy,iterMax);
        end

        function saveData(obj,step,uFun,dmgFun,qMax,rMax,uVal,fVal,energy,iterMax)
            s.uFun    = uFun;
            s.uVal    = uVal;
            s.fVal    = fVal;
            s.dmgFun  = dmgFun;
            s.qMax    = qMax;
            s.rMax    = rMax;
            s.energy  = energy;
            s.numIter = iterMax;
            obj.monitor.saveData(step,s)
        end

    end
end