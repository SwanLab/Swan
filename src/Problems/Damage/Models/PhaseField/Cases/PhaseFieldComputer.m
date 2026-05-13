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
            theta = obj.initialGuess.theta;
            cost = 0; tauArray = [];
            
            step = 0; iterMax.stag = 0;
            bc.phi   = obj.boundaryConditions.phi.updateStep(1); %Only done once to set initial damage field
            bcUfinal = obj.boundaryConditions.u.endVal;
            bcU      = obj.boundaryConditions.u.initVal;
            while(bcU <= bcUfinal) && (obj.stop.noFailure)
                notConverged = true; 
                maxIterReached = false; 
                convergenceTry = 1;
                while notConverged && (convergenceTry < 20)
                    [step,u,bc] = obj.updateBoundaryConditions(u,bc,step,iterMax.stag,maxIterReached);
                    bcU = obj.boundaryConditions.u.currentVal;
                    obj.monitor.printStep(bcU,bcUfinal,step,maxIterReached)
    
                    [u,theta,phi,F,cost,iterMax,maxIterReached] = obj.optimizer.compute(u,theta,phi,bc,cost);
                    convergenceTry = convergenceTry + 1;

                    if maxIterReached == true; notConverged = true;
                    else; notConverged = false; end
                end
                [Evec,totE,totF,uBC,angleVec,phiRel] = obj.postprocess(u,phi,theta,F,bc);
                obj.printAndSave(step,totF,uBC,u,theta,angleVec,phi,phiRel,Evec,totE,iterMax,cost,tauArray);
                obj.checkStopCondition(step,totF);

                % sig = obj.functional.computeStress(u,phi);
                % max(sig.evaluate([0;0]),[],'all')

                % eFun = obj.functional.computeEnergyFun(u);
                % gDeriv = obj.functional.computeDegradationDerivative(phi);
                % figure(2)
                % eFun.plot;
                % figure(3)
                % gDeriv.plot;
            end
            % if you want to retrieve the data at any point of the
            % simulation, just save the obj.monitor.data at any step!
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

            thetaPseudoFun.getDofCoord = obj.mesh.computeBaricenter';
            thetaPseudoFun.mesh.nnodes = size(thetaPseudoFun.getDofCoord,1);
            s.funs         = [{obj.initialGuess.phi.fun},{thetaPseudoFun}];
            s.legends      = [{["PhiMax","PhiRel"]}];
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

        function [step,u,bc] = updateBoundaryConditions(obj,u,bc,step,iterStag,maxIterReached)
            [bc.u,step] = obj.boundaryConditions.u.updateStep(step,iterStag,maxIterReached);
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

        function [E,totE,totF,uBC,angleVec,phiRel] = postprocess(obj,u,phi,theta,F,bc)
            fExt = bc.u.tractionFun;
            if ~isempty(bc.u.tractionFun)
                vals = bc.u.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.nnodes,u.mesh.ndim));
            end
            E    = obj.functional.computeEnergies(u,phi,fExt);
            totE = sum(E);
            [totF,uBC] = obj.computeTotalReaction(F,u);
            angleVec = obj.computeOrientationAsVector(theta,phi);
            phiRel   = obj.computeRelativeDamage(phi);
        end

        function [totReact,uBC] = computeTotalReaction(obj,F,u)
            DownSide = min(obj.mesh.coord(:,2));
            LeftSide = min(obj.mesh.coord(:,1));
            isInDown = abs(obj.mesh.coord(:,2)-DownSide)< 1e-12;
            isInLeft = abs(obj.mesh.coord(:,1)-LeftSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;

            if ismember(obj.boundaryConditions.u.type, ["ForceTractionY", "ForceTractionYClamped"])
                uBC = norm(mean(u.fValues(nodes(isInUp),2)));
                totReact = obj.boundaryConditions.u.currentVal;
            elseif ismember(obj.boundaryConditions.u.type, ["DisplacementTractionY","DisplacementTractionYClamped"]) 
                dofsYdown = (nodes(isInDown)-1)*u.ndimf + 2;
                totReact = abs(sum(F(dofsYdown)));
                uBC = obj.boundaryConditions.u.currentVal;
            end

            if ismember(obj.boundaryConditions.u.type, ["ForceTractionX","ForceTractionXClamped"])
                uBC = norm(mean(u.fValues(nodes(isInLeft),2)));
                totReact = obj.boundaryConditions.u.currentVal;
            elseif ismember(obj.boundaryConditions.u.type, ["DisplacementTractionX","DisplacementTractionXClamped"])
                dofsXleft = (nodes(isInLeft)-1)*u.ndimf + 1;
                totReact = abs(sum(F(dofsXleft)));
                uBC = obj.boundaryConditions.u.currentVal;
            end

            if ismember(obj.boundaryConditions.u.type, "DisplacementShear")
                dofsXdown = (nodes(isInDown)-1)*u.ndimf + 1;
                totReact = abs(sum(F(dofsXdown)));
                uBC = obj.boundaryConditions.u.currentVal;
            end

            if ismember(obj.boundaryConditions.u.type, "DisplacementPureShear")
                dofsCorner = (nodes(isInDown & isInLeft)-1)*u.ndimf + [1,2];
                totReact = norm(F(dofsCorner));
                uBC = obj.boundaryConditions.u.currentVal;
            end
        end

        function thetaVectorScaled = computeOrientationAsVector(~,theta,phi)
            thetaVal = squeeze(theta.evaluate([1/3;1/3])); % Computed at baricenter!
            phiVal = squeeze(phi.fun.evaluate([1/3;1/3]));
            thetaVector = [cos(thetaVal), sin(thetaVal)];
            thetaVectorScaled = phiVal.*thetaVector;
        end

        function phiRel = computeRelativeDamage(~,phi)
            maxPhi = max(phi.fun.fValues);
            AvgPhi = Integrator.compute(phi.fun,phi.mesh,2);
            phiRel = maxPhi/(AvgPhi + 1e-20);
        end

        function printAndSave(obj,step,totF,uBC,u,theta,angleVec,phi,phiRel,Evec,totE,iterMax,cost,tauArray)
            obj.monitor.updateAndRefresh(step,{[totF;uBC],[max(phi.fun.fValues);phiRel;uBC],...
                                               [phi.fun.fValues],[angleVec], ...
                                               [iterMax.stag],[],[totE;uBC],[]});
            obj.saveData(step,totF,uBC,u,theta,phi,Evec,iterMax,cost,tauArray);
        end

        function saveData(obj,step,totF,uVal,u,theta,phi,Evec,iterMax,cost,tauArray)
            s.force    = totF;
            s.bcVal    = uVal;
            s.u        = u;
            s.theta    = theta;
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
            
            if  (obj.stop.stepTrigger~=0) && (step==obj.stop.stepTrigger+10)
                obj.stop.noFailure = false;
            end
        end

    end
end
