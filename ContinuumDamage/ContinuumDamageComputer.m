classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material
        solverParams
    end

    properties (Access = private)
        tau = 0.5
        tolerance = 1e-8
        quadOrder

        H = 0.5
        r0 = 1/sqrt(210) %revisar com es calcula (depen de les bc)

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
                data.displacement.value(i)  = obj.boundaryConditions.bcValueSet(i);
                data.totalEnergy(i) = energy;
                damageDomainFun = obj.functional.computeDamage();
                damageFun = damageDomainFun.project('P1D',obj.mesh);
                data.damage.maxValue(i)  = max(damageFun.fValues);
                data.reaction(i)  = -obj.computeTotalReaction(LHS,uNewVec);
            end
            data.displacement.field = u;
            data.damage.field = damageFun;
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

        function totReact = computeTotalReaction(obj,LHS,u)
            F = LHS*u;
            isInDown =  (abs(obj.mesh.coord(:,2) - min(obj.mesh.coord(:,2)))< 1e-12);
            nodes = 1:obj.mesh.nnodes;
            totReact = sum(F(2*nodes(isInDown)));
        end

    end
end