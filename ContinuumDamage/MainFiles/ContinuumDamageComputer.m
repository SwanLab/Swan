classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        solverParams
        tolerance
        quadOrder
        internalDamageVariable
    end

    properties (Access = private)
        damageFunctional
        tau = 5e-4;
        limIter = 100;
    end

    methods (Access = public)
        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
        end

        function data = compute(obj)
            uFun = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.damageFunctional.setTestFunctions(uFun);
            data = {};

            nSteps = length(obj.boundaryConditions.bcValueSet);
            for i = 1:nSteps
                fprintf('Step: %d / %d \n',i,nSteps);
                bc = obj.updateBoundaryConditions(i);
                uFun.setFValues(obj.updateInitialDisplacement(bc,uFun));

                err = 1; iter = 0;
                while (err >= obj.tolerance && iter < obj.limIter)
                    [res,Ksec,uVec,uFun] = obj.solveU(uFun,bc);

                    err = norm(res);
                    fprintf('Error: %d \n',err);

                    iter = iter+1;
                end
                if (iter >= obj.limIter)
                    fprintf (2,'NOT CONVERGED FOR STEP %d\n',i);
                end

                obj.internalDamageVariable.updateRold();
                [data,dmgFun,rFun,qFun] = obj.getData(data,i,Ksec,uVec,uFun,bc);

            end
            data.displacement.field = uFun;
            data.damage.field = dmgFun;
            data.r.field = rFun;
            data.q.field = qFun;

        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.tolerance          = cParams.tol;
            obj.damageFunctional   = cParams.damageFunctional;
            obj.internalDamageVariable = cParams.internalDamageVariable;
        end

        function [bc] = updateBoundaryConditions (obj,i)
            bc = obj.boundaryConditions.nextStep(i);
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

        function [res,Ksec,uVec,u] = solveU(obj,u,bc)   
            tauEps = obj.damageFunctional.computeTauEpsilon(u);
            obj.internalDamageVariable.update(tauEps);
            r = obj.internalDamageVariable;
            [res]           = obj.damageFunctional.computeResidual(u,r,bc);
            [resDeriv,Ksec] = obj.damageFunctional.computeDerivativeResidual(u,r,bc);
            [uNew,uVec]     = obj.computeDisplacement(resDeriv,res,u,bc);
            u.setFValues(uNew);
            [res]  = obj.damageFunctional.computeResidual(u,bc);
        end

        function [uOut,uOutVec] = computeDisplacement(obj,LHS,RHS,uIn,bc)
            uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
            uOutVec = uInVec;

            uInFree = uInVec(bc.free_dofs);
            uOutFree = obj.updateWithNewton(LHS,RHS,uInFree);
            % uOutFree = obj.updateWithGradient(RHS,uInFree);
            uOutVec(bc.free_dofs) = uOutFree;
            uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';
        end

        % function xNew = updateWithGradient(obj,RHS,x)
        %     deltaX = -obj.tau.*RHS;
        %     xNew = x + deltaX;
        % end

        function xNew = updateWithNewton(~,LHS,RHS,x)
            deltaX = -LHS\RHS;
            xNew = x + deltaX;
        end

        function [data,dmgFun,rFun,qFun] = getData(obj,data,i,Ksec,uVec,uFun,bc)
            data.displacement.value(i)  = obj.boundaryConditions.bcValueSet(i);
            dmgDomainFun = obj.damageFunctional.getDamage();
            dmgFun = dmgDomainFun.project('P1D');
            data.damage.maxValue(i)  = max(dmgFun.fValues);
            data.damage.minValue(i)  = min(dmgFun.fValues);

            rDomainFun = obj.damageFunctional.getR();
            rFun = rDomainFun.project('P0');
            data.r.maxValue(i) = max(rFun.fValues);
            data.r.minValue(i) = min(rFun.fValues);

            qDomainFun = obj.damageFunctional.getQ();
            qFun = qDomainFun.project('P1D');
            data.q.maxValue(i) = max(qFun.fValues);
            data.q.minValue(i) = min(qFun.fValues);

            data.reaction(i)  = obj.computeTotalReaction(Ksec,uVec);
            [data.totalEnergy(i),data.damagedMaterial(i)] = obj.damageFunctional.computeEnergy(uFun,bc);
        end

        function totReact = computeTotalReaction(obj,LHS,u)
            F = LHS*u;
            isInDown =  (abs(obj.mesh.coord(:,2) - min(obj.mesh.coord(:,2)))< 1e-12);
            nodes = 1:obj.mesh.nnodes;
            totReact = -sum(F(2*nodes(isInDown)));
        end

    end
end