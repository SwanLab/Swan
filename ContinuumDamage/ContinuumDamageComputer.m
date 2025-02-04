classdef ContinuumDamageComputer < handle

    properties (Access = private)
        tolerance

        mesh
        boundaryConditions
        material
        H
        solverParams
        quadOrder 
 
    end

    properties (Access = private)
        tau = 5e-5;
        r0

        elasticFun
        externalWorkFun
        elasticity
    end

    methods (Access = public)
        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
            obj.defineRfunction(cParams)
            obj.createFunctionals()
        end

        function data = compute(obj)
            uFun = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.elasticity.setTestFunctions(uFun);
            
            for i = 1:obj.boundaryConditions.valueSetLenght
                fprintf('Step: %d ',i);fprintf('/ %d \n',obj.boundaryConditions.valueSetLenght);
                bc = obj.updateBoundaryConditions(i);
                uFun.setFValues(obj.updateInitialDisplacement(bc,uFun));

                resErr = 1;
                while (resErr >= obj.tolerance)
                    obj.elasticity.computeDamageEvolutionParam(uFun);
                    [res]  = obj.elasticity.computeResidual(uFun,bc);
                    [K,resDeriv] = obj.elasticity.computeDerivativeResidual(uFun,bc);
                    [uVal,uVec] = obj.computeDisplacement(resDeriv,res,uFun,bc);
                    uFun.setFValues(uVal);

                    resErr = norm(res);
                    fprintf('Error: %d \n',resErr);
                end
                obj.elasticity.setROld();

                data.displacement.value(i)  = obj.boundaryConditions.bcValueSet(i);
                dmgDomainFun = obj.elasticity.getDamage();
                dmgFun = dmgDomainFun.project('P1D');
                data.damage.maxValue(i)  = max(dmgFun.fValues);
                data.damage.minValue(i)  = min(dmgFun.fValues);
                data.reaction(i)  = -obj.computeTotalReaction(K,uVec);
            end
            data.displacement.field = uFun;
            data.damage.field = dmgFun;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material  = cParams.material;
            obj.solverParams = cParams.solver; 
            obj.H = cParams.H;
            obj.tolerance = cParams.tol;
            obj.quadOrder = 2; %maybe try to set it from Main?
        end

        function defineRfunction(obj,cParams)
            obj.r0 = LagrangianFunction.create(obj.mesh,1,'P0');
            fV = cParams.r0*ones(size(obj.r0.fValues));
            obj.r0.setFValues(fV);
        end
        
        function bc = updateBoundaryConditions (obj,i)
            bc = obj.boundaryConditions.nextStep(i);
        end

        function createFunctionals(obj)
            s.mesh     = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.material = obj.material;
            s.H = obj.H;
            s.r0 = obj.r0;
            s.quadOrder = obj.quadOrder;
            obj.elasticity = ShFunc_ContinuumDamage(s);
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

        function [uOut,uOutVec] = computeDisplacement(obj,LHS,RHS,uIn,bc)            
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