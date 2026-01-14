classdef PhaseFieldComputer < handle

    properties (Access = private)
        mesh
        initialGuess
        boundaryConditions
        functional
    end

    properties (Access = private)
        monitor
        optimizer
        stop
    end

    methods (Access = public)

        function obj = PhaseFieldComputer(cParams)
            obj.init(cParams)
            obj.setMonitoring(cParams)
            obj.setOptimizer(cParams)
            obj.setStopConditions()
        end

        function outputData = compute(obj)
            u   = obj.initialGuess.u;
            phi = obj.initialGuess.phi;
            cost = 0; tauArray = [];

            step = 1;
            maxSteps = length(obj.boundaryConditions.u.bcValues);
            bc.phi = obj.boundaryConditions.phi.nextStep();
            while(step<=maxSteps) && (obj.stop.noFailure)
                obj.monitor.printStep(step,maxSteps)
                [u,bc] = obj.updateBoundaryConditions(u,bc);
                [u,phi,F,cost,iterMax] = obj.optimizer.compute(u,phi,bc,cost);
                [Evec,totE,totF,uBC] = obj.postprocess(step,u,phi,F,bc);
                obj.printAndSave(step,totF,uBC,u,phi,Evec,totE,iterMax,cost,tauArray);
                obj.checkStopCondition(step,totF);
                step = step + 1;

                sig = obj.functional.computeStress(u,phi);
                max(sig.evaluate([0;0]),[],'all')
            end
            outputData = obj.monitor.data;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.initialGuess       = cParams.initialGuess;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.functional         = cParams.functional;
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay = cParams.monitoring.set;
            s.shallPrint   = cParams.monitoring.print;
            s.type         = cParams.monitoring.type;
            s.fun          = obj.initialGuess.phi.fun;
            obj.monitor = PhaseFieldMonitoring(s);
        end

        function setOptimizer(obj,cParams)
            s.functional  = obj.functional;
            s.initPhi     = obj.initialGuess.phi;
            s.tolerance   = cParams.tolerance;
            s.maxIter     = cParams.maxIter;
            s.solver      = cParams.solver;
            s.monitor     = obj.monitor;
            obj.optimizer = OptimizerPhaseField(s);
        end

        function setStopConditions(obj)
            maxSteps = length(obj.boundaryConditions.u.bcValues);
            obj.stop.noFailure   = true;
            obj.stop.triggered   = false;
            obj.stop.maxF        = 0;
            obj.stop.stepTrigger = maxSteps;
        end

        function [u,bc] = updateBoundaryConditions(obj,u,bc)
            bc.u = obj.boundaryConditions.u.nextStep();
            u.setFValues(obj.updateInitialDisplacement(u,bc));
        end

        function u = updateInitialDisplacement(~,uOld,bc)
            restrictedDofs = bc.u.dirichlet_dofs;
            if isempty(restrictedDofs)
                u = uOld;
            else
                dirich = bc.u.dirichletFun;
                uVec = reshape(uOld.fValues',[uOld.nDofs 1]);
                dirichVec = reshape(dirich.fValues',[dirich.nDofs 1]);

                uVec(restrictedDofs) = dirichVec(restrictedDofs);
                u = reshape(uVec,[flip(size(uOld.fValues))])';
            end
        end

        function [E,totE,totF,uBC] = postprocess(obj,step,u,phi,F,bc)
            fExt = bc.u.tractionFun;
            if ~isempty(bc.u.tractionFun)
                vals = bc.u.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.nnodes,u.mesh.ndim));
            end
            E    = obj.functional.computeEnergies(u,phi,fExt);
            totE = sum(E);
            [totF,uBC] = obj.computeTotalReaction(step,F,u);
        end

        function [totReact,uBC] = computeTotalReaction(obj,step,F,u)
            LeftSide = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-LeftSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            if ismember(obj.boundaryConditions.u.type, ["ForceTractionY", "ForceTractionYClamped"])
                uBC = norm(mean(u.fValues(nodes(isInUp),2)));
                totReact = obj.boundaryConditions.u.bcValues(step);
            elseif ismember(obj.boundaryConditions.u.type, ["DisplacementTractionY","DisplacementTractionYClamped"]) 
                totReact = norm(sum(F(2*nodes(isInUp))));
                uBC = obj.boundaryConditions.u.bcValues(step);
            end

            LeftSide = min(obj.mesh.coord(:,1));
            isInLeft = abs(obj.mesh.coord(:,1)-LeftSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            if ismember(obj.boundaryConditions.u.type, ["ForceTractionX","ForceTractionXClamped"])
                uBC = norm(mean(u.fValues(nodes(isInLeft),2)));
                totReact = obj.boundaryConditions.u.bcValues(step);
            elseif ismember(obj.boundaryConditions.u.type, ["DisplacementTractionX","DisplacementTractionXClamped"])
                dofsXleft = (nodes(isInLeft)-1)*u.ndimf + 1;
                totReact = abs(sum(F(dofsXleft)))/10;
                uBC = obj.boundaryConditions.u.bcValues(step);
            end
            

            % DownSide  = min(obj.mesh.coord(:,2));
            % isInDown = abs(obj.mesh.coord(:,2)-DownSide)< 1e-12;
            % nodes = 1:obj.mesh.nnodes;
            % totReact = -sum(F(2*nodes(isInDown)));

            % isInTip = (abs(obj.mesh.coord(:,2)-(max(obj.mesh.coord(:,2))+min(obj.mesh.coord(:,2)))/2) < 1e-12) & (abs(obj.mesh.coord(:,1)-max(obj.mesh.coord(:,1))) < 30);
            % nodes = 1:obj.mesh.nnodes;
            % totReact = sum(F(2*nodes(isInTip)));
        end

        function printAndSave(obj,step,totF,uBC,u,phi,Evec,totE,iterMax,cost,tauArray)
            obj.monitor.updateAndRefresh(step,{[totF;uBC],[max(phi.fun.fValues);uBC],...
                                               [phi.fun.fValues],[iterMax.stag],[],...
                                               [totE;uBC],[]});
            obj.saveData(step,totF,uBC,u,phi,Evec,iterMax,cost,tauArray);
        end

        function saveData(obj,step,totF,uVal,u,phi,Evec,iterMax,cost,tauArray)
            s.force    = totF;
            s.bcVal    = uVal;
            s.u        = u;
            s.phi      = phi;
            s.energy   = Evec;
            s.numIter  = iterMax;
            s.cost     = cost;
            s.tauArray = tauArray;
            obj.monitor.saveData(step,s);
        end

        function checkStopCondition(obj,step,totF)
            if totF > obj.stop.maxF
                obj.stop.maxF = totF;
            elseif step>5 && totF<0.01*obj.stop.maxF && ~obj.stop.triggered
                obj.stop.stepTrigger = step;
                obj.stop.triggered = true;
            end

            if step==obj.stop.stepTrigger+10
                obj.stop.noFailure = false;
            end
        end

    end
end
