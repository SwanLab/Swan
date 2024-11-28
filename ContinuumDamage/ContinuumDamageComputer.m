classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material
        solverParams
    end

    properties (Access = private)
        tau = 0.1
        tolerance = 1e-5
        quadOrder

        H = 0.1
        r0 = (4.0e-1)/sqrt(3e4)%revisar com es calcula (depen de les bc)

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
            bc = obj.boundaryConditions;
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            u.fValues = obj.updateInitialDisplacement(bc,u);
                        
            errorE = 1;
            fExt = obj.boundaryConditions.pointloadFun;
            EnergyOld = -1;

            while (abs(errorE) >= obj.tolerance)
                LHS = obj.functional.computeHessian(obj.quadOrder,u);
                RHS = obj.functional.computeJacobian(obj.quadOrder,u,bc);
                [uNew,uNewVec] = obj.computeU(LHS,RHS,u,bc);
                EnergyNew = obj.functional.computeTotalEnergyDamage(obj.quadOrder,u,fExt);

                errorE = max(max(EnergyNew-EnergyOld));
                fprintf('Error: %d ',errorE);
                fprintf('Cost: %d \n',EnergyNew);

                EnergyOld = EnergyNew;
                u.fValues = uNew;
                obj.functional.updateInternalVariableR(u);
            end
            data.displacement = u;
            data.displacement.value = obj.boundaryConditions.bcVals;
            data.TotalEenrgy = EnergyNew;
            damageDomainFun = obj.functional.computeDamage();
            damageFun = damageDomainFun.project('P1D',obj.mesh);
            data.damage.fun = damageFun;
            data.damage.maxValue = max(damageFun.fValues);
            data.reactions = obj.computeReactions(u,uNewVec,LHS);
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

            obj.elasticFun = shFunc_ElasticDamage(s);
            obj.externalWorkFun = shFunc_ExternalWork2(s);
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