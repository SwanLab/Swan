classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material
        solverParams
    end

    properties (Access = private)
        tau = 50
        tolerance = 1e-10
        quadOrder

        H = 0.01
        r0 = 0.1

        ElasticFun
        ExternalWorkFun
        TotalEnergyFun
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
            
            rNew = ConstantFunction.create(obj.r0,obj.mesh);
            
            errorU = 1;
            while (errorU >= obj.tolerance)
                LHS = obj.computeLHS(u,rNew);
                RHS = obj.computeRHS(u,rNew);
                uNew = obj.computeU(LHS,RHS,u,bc);
                errorU = max(max(abs(u.fValues-uNew)));
                
                u.fValues = uNew;
                rOld = rNew;
                rNew = obj.ElasticFun.newState(rOld,u);

                fprintf('Error: %d \n',errorU);
            
            end
            data.displacement = u;
            fExt = obj.boundaryConditions.pointloadFun;
            %data.TotalEenrgy = obj.TotalEnergyFun.computeTotalEnergy(obj.quadOrder,u,fExt);
            data.TotalEenrgy = obj.TotalEnergyFun.computeTotalEnergyDamage(obj.quadOrder,u,rNew,fExt);
            data.damage = obj.ElasticFun.computeDamage(rNew);
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

            obj.ElasticFun = shFunc_ElasticDamage(s);
            %obj.ElasticFun = shFunc_Elastic(s);

            obj.ExternalWorkFun = shFunc_ExternalWork2(s);
            obj.TotalEnergyFun = shFunc_TotalEnergy(s)
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
            K = obj.ElasticFun.computeHessian(obj.quadOrder,u,r);
        end

        function F = computeRHS(obj,u,r)
            fExt = obj.boundaryConditions.pointloadFun;
            Fext = obj.ExternalWorkFun.computeGradient(u,fExt,obj.quadOrder);
            Fint = obj.ElasticFun.computeJacobian(obj.quadOrder,u,r);
            F = Fint - Fext;
        end

        function uOut = computeU(obj,LHS,RHS,uIn,bc)            
            RHS = RHS(bc.free_dofs);
            LHS = LHS(bc.free_dofs,bc.free_dofs);

            uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
            uOutVec = uInVec;

            uInFree = uInVec(bc.free_dofs);
            %uOutFree = obj.updateWithNewton(LHS,RHS,uInFree);
            uOutFree = obj.updateWithGradient(RHS,uInFree);
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