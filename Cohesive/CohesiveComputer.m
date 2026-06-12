classdef CohesiveComputer < handle
    
    properties (Access = public) 
        mesh
        boundaryConditions       
        functional
        tolerance
        maxIter
        data
    end
    
    properties (Access = private)
        optimizer
        cohesiveMesh
        updater
        jump
        tractionSeparation
        monitor
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveComputer(cParams)
            obj.init(cParams)
            obj.setMonitoring(cParams);
            obj.setOptimizer(cParams)
        end

        function compute(obj)
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            cost = 0;

            nSteps = length(obj.boundaryConditions.bcValues);
            for iStep = 1:nSteps
                
                [u,bc] = obj.preprocess(iStep,nSteps,u);
                [u,F,cost,iterMax] = obj.updater.update(u,bc,cost);
                obj.postprocess(iStep,u,F,cost,iterMax)
            end
            obj.data = obj.data;
            outputData = obj.monitor.data;
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh               = cParams.cohesiveMesh.fullMesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.functional         = cParams.functional;      
            obj.tractionSeparation = cParams.tractionLaw;
            obj.cohesiveMesh       = cParams.cohesiveMesh;
        end
        
        function setOptimizer(obj,cParams)
            s.functional = obj.functional;
            s.tolerance  = cParams.tolerance;
            s.maxIter    = cParams.maxIter;
            s.solverType = cParams.solverType;
            s.monitor    = 1;
            obj.updater = DisplacementUpdater(s);
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay   = cParams.monitoring.set;
            s.shallPrintInfo = cParams.monitoring.print;
            s.fun            = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.monitor = CohesiveMonitoring(s);
        end

        function [u,bc] = preprocess(obj,iStep,nSteps,u)
            obj.monitor.printStep(iStep,nSteps)
            bc = obj.boundaryConditions.nextStep();
            u  = obj.computeInitialDisplacement(u,bc);
        end

        function postprocess(obj,iStep,uFun,F,cost,iterMax)
            [fVal,uVal] = obj.computeTotalReaction(iStep,F,uFun);
            dFun = obj.computeDamageField(uFun);
            obj.printAndSave(iStep,uFun,dFun,uVal,fVal,cost(end),iterMax);
            % obj.data = [obj.data; fVal,dFun,uFun,uVal];
            obj.data = [obj.data; fVal,uVal, dFun.fValues];
        end

        function [totReact,uBC] = computeTotalReaction(obj,step,F,u)
            DownSide = min(obj.mesh.coord(:,2));
            isInDown = abs(obj.mesh.coord(:,2)-DownSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            if  ismember(obj.boundaryConditions.type, ["DisplacementTractionY","DisplacementTractionYClamped"]) 
                dofsYdown = (nodes(isInDown)-1)*u.ndimf + 2;
                totReact = abs(sum(F(dofsYdown)));
                uBC = obj.boundaryConditions.bcValues(step);
            end
           isRight   = abs(obj.mesh.coord(:,1) - max(obj.mesh.coord(:,1))) < 1e-12;
           MiddleY   = (max(obj.mesh.coord(:,2)) + min(obj.mesh.coord(:,2)))/2;
           isHalfTop = (obj.mesh.coord(:,2) - MiddleY) > 1e-12;
           isRigthHalfTop = isHalfTop & isRight;
            if ismember(obj.boundaryConditions.type, "DoubleCantileverBeam")
                dofYHalfRight = (nodes(isRigthHalfTop)-1)*u.ndimf + 2;
                totReact = abs(sum(F(dofYHalfRight)));
                uBC = obj.boundaryConditions.bcValues(step);
            end
            if ismember(obj.boundaryConditions.type, "EndNotchedFlex")
               isInBottomLeft =  (abs(obj.mesh.coord(:,1) - min(obj.mesh.coord(:,1)))< 1e-12) & (abs(obj.mesh.coord(:,2) - min(obj.mesh.coord(:,2)))< 1e-12);
               isInBottomRight = (abs(obj.mesh.coord(:,1) - min(obj.mesh.coord(:,1)))< 1e-12) & (abs(obj.mesh.coord(:,2) - min(obj.mesh.coord(:,2)))< 1e-12);
               nodes = find(isInBottomRight | isInBottomLeft);
               dofsY= (nodes-1)*u.ndimf + 2;
               totReact = abs(sum(F(dofsY)));
               uBC = obj.boundaryConditions.bcValues(step);
            end
        end

        function u = computeInitialDisplacement(obj,u,bc)
            uVec = reshape(u.fValues',[u.nDofs 1]);
            uVec(bc.dirichlet_dofs) = bc.dirichlet_vals;
            ufV = reshape(uVec,[flip(size(u.fValues))])';
            u.setFValues(ufV);
        end

        function dFun = computeDamageField(obj,u)
            s.cohesiveMesh = obj.cohesiveMesh;
            s.ndimf = obj.cohesiveMesh.fullMesh.ndim;
            s.uFun  = LagrangianFunction.create(obj.cohesiveMesh.fullMesh,s.ndimf,'P1');
            obj.jump = Jump(s);
            obj.jump.updateJumpValues(u);
            dValues = obj.tractionSeparation.getDamageValues(obj.jump.fun);
            dValues = dValues.evaluate([-1,1]);

            temp1 = reshape(dValues,1,[]);
            temp2 = temp1(2:end-1);
            temp3 = reshape(temp2,2,size(temp2,2)/2);
            temp4 = sum(temp3,1)*0.5;
            dValues = [temp1(1),temp4,temp1(end)];

            dFun = LagrangianFunction.create(obj.cohesiveMesh.mesh,size(dValues,1),'P1');
            dFun.setFValues(dValues);
        end

        function printAndSave(obj,iStep,uFun,dmgFun,uVal,fVal,energy,iterMax)
            dmgMax = max(dmgFun.fValues); 
            obj.monitor.updateAndRefresh(iStep,{[fVal;uVal],[dmgMax;uVal],...
                [energy],[iterMax]});
            obj.saveData(iStep,uFun,dmgFun,uVal,fVal,energy,iterMax);
        end


        function saveData(obj,step,uFun,dmgFun,uVal,fVal,energy,iterMax)
            s.uFun    = uFun;
            s.uVal    = uVal;
            s.fVal    = fVal;
            s.dmgFun  = dmgFun;
            s.energy  = energy;
            s.numIter = iterMax;
            obj.monitor.saveData(step,s)
        end

    end
    
end