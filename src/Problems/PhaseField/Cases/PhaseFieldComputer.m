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
        end

        function outputData = compute(obj)
            u   = obj.initialGuess.u;
            phi = obj.initialGuess.phi;
            
            uOld = u;
            cost = 0;
            tauArray = [];

            maxSteps = length(obj.boundaryConditions.bcValues);
            obj.stop.noFullyBroken = true; obj.stop.triggered = false;
            obj.stop.maxF = 0;             obj.stop.lastStep = maxSteps;
            i = 1;
            while(i<=maxSteps) && (obj.stop.noFullyBroken)
                obj.monitor.printStep(i,maxSteps)
                bc = obj.boundaryConditions.nextStep();
                u.setFValues(obj.updateInitialDisplacement(bc,uOld));
                [u,phi,F,cost,iterMax] = obj.optimizer.compute(u,phi,bc,cost);
                uOld = u;


                % if obj.solverType == "Gradient"
                %     % ADAPTIVE LINE-SEARCH GRADIENT
                %     phiNew = obj.updateWithGradient(RHS, phi.fValues,tau);
                %     phiProposed = phi.copy();
                %     phiProposed.setFValues(obj.projectInLowerAndUpperBound(phiNew,phiOld.fValues,1));
                %     [ePhi, costPhi] = obj.computeErrorCostFunctional(u,phiProposed,bc,costOldPhi);
                %
                %     obj.monitor.printCost('iterPhi',iterPhi,costPhi,ePhi);
                %     iterPhi = iterPhi + 1;
                %
                %     if ePhi > 0
                %         tau = tau/2; %% Estem el mateix RHS i LHS sempre que entrem aqui
                %     else
                %         phiIter = phiIter + 1;
                %         obj.monitor.update(phiIter,{[],[],[],[],[],[],[tau]})
                %         tauArray(end+1) = tau;
                %         if tau<=1e10
                %             tau = 10*tau;
                %         end
                %         phi = phiProposed;
                %         costOldPhi = costPhi;
                %         costFun(end+1) = costPhi;
                %         iter = iter+1;
                %         obj.monitor.update(iter,{[],[],[],[],[costFun(end)],[],[]});
                %     end

                fExt = bc.pointloadFun;
                Evec = obj.functional.computeEnergiesFunctional(u,phi,fExt);
                totE = sum(Evec);
                totF = obj.computeTotalReaction(F);

                uVal = obj.boundaryConditions.bcValues(i);
                obj.monitor.update(i,{[totF;uVal],[max(phi.fValues);uVal],[phi.fValues],...
                                      [iterMax.stag],[],[totE;uVal],[]});
                obj.monitor.refresh();

                s.force    = totF;
                s.bcVal    = obj.boundaryConditions.bcValues(i);
                s.u        = u;
                s.phi      = phi;
                s.energy   = Evec;
                s.numIter  = iterMax;
                s.cost     = cost;
                s.tauArray = tauArray;
                obj.monitor.saveData(i,s);


                obj.checkStopCode(i,totF);
                i = i + 1;
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
            obj.setMonitoring(cParams)
            obj.setOptimizer(cParams)
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay = cParams.monitoring.set;
            s.shallPrint   = cParams.monitoring.print;
            s.type         = cParams.monitoring.type;
            s.fun          = obj.initialGuess.phi;
            obj.monitor = PhaseFieldMonitoring(s);
        end

        function setOptimizer(obj,cParams)
            s.functional  = obj.functional;
            s.initPhi     = obj.initialGuess.phi;
            s.tolerance   = cParams.tolerance;
            s.maxIter     = cParams.maxIter;
            s.solverType  = cParams.solverType;
            s.monitor     = obj.monitor;
            obj.optimizer = OptimizerPhaseField(s);
        end

        function u = updateInitialDisplacement(obj,bc,uOld)
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

        function checkStopCode(obj,i,totF)
            if totF > obj.stop.maxF
                obj.stop.maxF = totF;
            elseif i>5 && totF<0.01*obj.stop.maxF && ~obj.stop.triggered
                obj.stop.lastStep = i;
                obj.stop.triggered = true;
            end

            if i==obj.stop.lastStep+10
                obj.stop.noFullyBroken = false;
            end
        end

        function totReact = computeTotalReaction(obj,F)
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            ReactX = sum(F(2*nodes(isInUp)-1));
            ReactY = sum(F(2*nodes(isInUp)));
            totReact = sqrt(ReactX^2+ReactY^2);

            % DownSide  = min(obj.mesh.coord(:,2));
            % isInDown = abs(obj.mesh.coord(:,2)-DownSide)< 1e-12;
            % nodes = 1:obj.mesh.nnodes;
            % totReact = -sum(F(2*nodes(isInDown)));

            % isInTip = (abs(obj.mesh.coord(:,2)-(max(obj.mesh.coord(:,2))+min(obj.mesh.coord(:,2)))/2) < 1e-12) & (abs(obj.mesh.coord(:,1)-max(obj.mesh.coord(:,1))) < 30);
            % nodes = 1:obj.mesh.nnodes;
            % totReact = sum(F(2*nodes(isInTip)));
        end

    end
end
