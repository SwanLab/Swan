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

        Functional
    end

    methods (Access = public)

        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
            obj.createFunctionals()
        end

        function data = compute(obj)
            bc = obj.boundaryConditions.bc;
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            u.fValues = obj.updateInitialDisplacement(bc,u);

            rNew = ConstantFunction.create(obj.r0,obj.mesh);

            errorE = 1;
            fExt = obj.boundaryConditions.bc.pointloadFun;

            EnergyOld = 1;
            
            

            for i = 1:1:obj.boundaryConditions.ValueSetLenght
                fprintf('Step: %d ',i);fprintf('/ %d \n',obj.boundaryConditions.ValueSetLenght);
                
                
                obj.boundaryConditions.nextStep(i);
                bc = obj.boundaryConditions.bc;
                
            

                while (abs(errorE) >= obj.tolerance)
                    LHS = obj.computeLHS(u,rNew);
                    RHS = obj.computeRHS(u,rNew);
                    [uNew,uNewVec] = obj.computeU(LHS,RHS,u,bc);

                    EnergyNew = obj.Functional.computeTotalEnergyDamage(obj.quadOrder,u,rNew,fExt);

                    errorE = max(max(EnergyNew-EnergyOld));

                    EnergyOld = EnergyNew;
                    u.fValues = uNew;
                    rOld = rNew;

                    rNew = obj.Functional.newState(rOld,u);
                    fprintf('Error: %d ',errorE);
                    fprintf('Cost: %d \n',EnergyNew);

                end
                errorE = 1;
                reactions(i) = obj.computeReactions (u,uNewVec,LHS);



            end
            data.displacement = u;
            %fInt = obj.Functional.computeJacobian(obj.quadOrder,u,rNew);
            %fExt = obj.boundaryConditions.pointloadFun;
            %data.TotalEenrgy = obj.TotalEnergyFun.computeTotalEnergy(obj.quadOrder,u,fExt);
            data.TotalEenrgy = EnergyNew;
            data.damage = obj.Functional.computeDamage(rNew);
            data.reactions = reactions;

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

            % obj.ElasticFun = shFunc_ElasticDamage(s);
            % obj.ElasticFun = shFunc_Elastic(s);
            %
            % obj.ExternalWorkFun = shFunc_ExternalWork2(s);
            obj.Functional = shFunc_ContinuumDamage(s);
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

        function K = computeLHS(obj,u,r)
            K = obj.Functional.computeHessian(obj.quadOrder,u,r);
        end

        function F = computeRHS(obj,u,r)
            fExt = obj.boundaryConditions.bc.pointloadFun;
            Fext = obj.Functional.computeGradient(u,fExt,obj.quadOrder);
            Fint = obj.Functional.computeJacobian(obj.quadOrder,u,r);
            F = Fint - Fext;
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

    end
end