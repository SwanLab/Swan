classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material
        solverParams
    end

    properties (Access = private)
        tau = 0.5
        tolerance = 1
        quadOrder

        H = 0.001
        r0 = 1e-3/sqrt(210) %revisar com es calcula (depen de les bc)

        elasticFun
        externalWorkFun
        functional
    end

    methods (Access = public)

        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
            obj.createFunctionals()
        end

        function data = compute(obj)
            u = LagrangianFunction.create(obj.mesh,2,'P1');

            for i = 1:obj.boundaryConditions.ValueSetLenght
                fprintf('Step: %d ',i);fprintf('/ %d \n',obj.boundaryConditions.ValueSetLenght);

                bc = obj.boundaryConditions.nextStep(i);
                fExt = bc.pointloadFun;
                u.fValues = obj.updateInitialDisplacement(bc,u);

                error = 1; energyOld = 1;
                while (abs(error) >= obj.tolerance)
                    LHS = obj.functional.computeHessian(obj.quadOrder,u);
                    RHS = obj.functional.computeJacobian(obj.quadOrder,u,bc);
                    [uNew,uNewVec] = obj.computeU(LHS,RHS,u,bc);
                    energy = obj.functional.computeTotalEnergyDamage(obj.quadOrder,u,fExt);

                    error = max(max(energy-energyOld));
                    fprintf('Error: %d ',error);
                    fprintf('Cost: %d \n',energy);

                    energyOld = energy;
                    u.fValues = uNew;
                    obj.functional.updateInternalVariableR(u);
                end
                data.displacement.value = obj.boundaryConditions.bcVals;
                data.totalEnergy = energy;
                damageDomainFun = obj.functional.computeDamage();
                damageFun = damageDomainFun.project('P1D',obj.mesh);
                data.damage.maxValue = max(damageFun.fValues);
                data.reactions = obj.computeReactions(u,uNewVec,LHS);
            end
            data.displacement.field = u;
            data.damage.fun = damageFun;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = 2;
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material  = cParams.material;
            obj.solverParams = cParams.solver;
        end

        function createFunctionals(obj)
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.H = obj.H;
            s.r0 = obj.r0;

            obj.functional = shFunc_ContinuumDamage(s);
        end

        function Reac = computeReactions (obj, uLagrangian, u, LHS)

            R = LHS*u;
            R(obj.boundaryConditions.bc.free_dofs) = 0;

            Rout = reshape(R,[flip(size(uLagrangian.fValues))])';
            Reac = LagrangianFunction.create(obj.mesh,2,'P1');

            Reac.fValues = Rout;

            isInDown =  (abs(obj.mesh.coord(:,2) - min(obj.mesh.coord(:,2)))< 1e-12);

            Reac.fValues(:,2) = Reac.fValues(:,2).*isInDown;
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

        function [uOut,uOutVec] = computeU(obj,LHS,RHS,uIn,bc)            
            RHS = RHS(bc.free_dofs);
            LHS = LHS(bc.free_dofs,bc.free_dofs);

            uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
            uOutVec = uInVec;

            uInFree = uInVec(bc.free_dofs);
            uOutFree = obj.updateWithNewton(LHS,RHS,uInFree);
            %uOutFree = obj.updateWithGradient(RHS,uInFree);
            uOutVec(bc.free_dofs) = uOutFree;
            uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';
        end

        function xNew = updateWithGradient(obj,RHS,x)
            deltaX = -obj.tau.*RHS;
            xNew = x + deltaX;
        end

        function xNew = updateWithNewton(~,LHS,RHS,x)
            deltaX = -LHS\RHS;
            xNew = x + deltaX;
        end

        function Reac = computeReactions (obj, uLagrangian, u, LHS)
            R = LHS*u;
            R(obj.boundaryConditions.free_dofs) = 0;
            
            Rout = reshape(R,[flip(size(uLagrangian.fValues))])';
            Reac = LagrangianFunction.create(obj.mesh,2,'P1');
            
            Reac.fValues = Rout;

            isInDown =  (abs(obj.mesh.coord(:,2) - min(obj.mesh.coord(:,2)))< 1e-12);
            Reac.fValues = Rout.*isInDown;
        end

    end
end