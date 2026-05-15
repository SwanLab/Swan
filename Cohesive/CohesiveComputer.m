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
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveComputer(cParams)
            obj.init(cParams)
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
            % outputData = obj.monitor.data;
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

        function [u,bc] = preprocess(obj,iStep,nSteps,u)
            % obj.monitor.printStep(iStep,nSteps)
            bc = obj.boundaryConditions.nextStep();
            u  = obj.computeInitialDisplacement(u,bc);
        end

        function postprocess(obj,iStep,uFun,F,cost,iterMax)
            [fVal,uVal] = obj.computeTotalReaction(iStep,F,uFun);
            % obj.printAndSave(iStep,uFun,dmgFun0,dmgFun1,qFun,rFun,uVal,fVal,cost(end),iterMax);
            % guardar
            dFun = obj.computeDamageField(uFun);
            obj.data = [obj.data; fVal,dFun,uFun,uVal];
            obj.data = [obj.data; uFun,uVal];
        end

        function [totReact,uBC] = computeTotalReaction(obj,step,F,u)
            F = reshape(F.fValues',[F.nDofs 1]);
            DownSide = min(obj.mesh.coord(:,2));
            isInDown = abs(obj.mesh.coord(:,2)-DownSide)< 1e-12;
            isInUp   = not(isInDown);
            nodes = 1:obj.mesh.nnodes;
            if ismember(obj.boundaryConditions.type, ["ForceTractionY", "ForceTractionYClamped"])
                uBC = norm(mean(u.fValues(nodes(isInUp),2)));
                totReact = obj.boundaryConditions.bcValues(step);
            elseif ismember(obj.boundaryConditions.type, ["DisplacementTractionY","DisplacementTractionYClamped"]) 
                dofsYdown = (nodes(isInDown)-1)*u.ndimf + 2;
                totReact = abs(sum(F(dofsYdown)));
                uBC = obj.boundaryConditions.bcValues(step);
            elseif ismember(obj.boundaryConditions.type, ["DisplacementMixed"])
                dofsXdown = (nodes(isInDown)-1)*u.ndimf + 1;
                dofsYdown = (nodes(isInDown)-1)*u.ndimf + 2;
                normF = sqrt(F(dofsXdown).^2+ F(dofsYdown).^2);
                totReact = abs(sum(normF));
                uBC = obj.boundaryConditions.bcValues(step);
            end

            LeftSide = min(obj.mesh.coord(:,1));
            isInLeft = abs(obj.mesh.coord(:,1)-LeftSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            if ismember(obj.boundaryConditions.type, ["ForceTractionX","ForceTractionXClamped"])
                uBC = norm(mean(u.fValues(nodes(isInLeft),2)));
                totReact = obj.boundaryConditions.bcValues(step);
            elseif ismember(obj.boundaryConditions.type, ["DisplacementTractionX","DisplacementTractionXClamped"])
                dofsXleft = (nodes(isInLeft)-1)*u.ndimf + 1;
                totReact = abs(sum(F(dofsXleft)));
                uBC = obj.boundaryConditions.bcValues(step);
            end

            if ismember(obj.boundaryConditions.type, "DisplacementShear")
                dofsXdown = (nodes(isInDown)-1)*u.ndimf + 1;
                totReact = abs(sum(F(dofsXdown)));
                uBC = obj.boundaryConditions.bcValues(step);
            end

            if ismember(obj.boundaryConditions.type, "DoubleCantileverBeam")
                idxNode = find(obj.mesh.coord(:,1) == max(obj.mesh.coord(:,1)) & (obj.mesh.coord(:,2) == max(obj.mesh.coord(:,2))));
                dofYTopRight = (idxNode-1)*u.ndimf + 1;
                totReact = abs(F(dofYTopRight));
                uBC = obj.boundaryConditions.bcValues(step);
            end

            if ismember(obj.boundaryConditions.type, "EndNotchedFlexural")
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

        function d = computeDamageField(obj,u)
            s.cohesiveMesh = obj.cohesiveMesh;
            s.ndimf = obj.cohesiveMesh.fullMesh.ndim;
            s.uFun  = LagrangianFunction.create(obj.cohesiveMesh.fullMesh,s.ndimf,'P1');
            obj.jump = Jump(s);
            obj.jump.updateJumpValues(u);
            d = obj.tractionSeparation.law.computeDamage(obj.jump.fun);
        end
        
    end
    
end