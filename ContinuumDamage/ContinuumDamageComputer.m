classdef ContinuumDamageComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        material
        solverParams
    end

    properties (Access = private)
        tau = 5e-5
        tolerance = 1e-8
        quadOrder

        H
        r0

        elasticFun
        externalWorkFun
        elasticity
    end

    methods (Access = public)

        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
            obj.createFunctionals()
        end

        function data = compute(obj)
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.elasticity.createTest(u);
            
            for i = 1:obj.boundaryConditions.ValueSetLenght
                fprintf('Step: %d ',i);fprintf('/ %d \n',obj.boundaryConditions.ValueSetLenght);

                bc = obj.updateBoundaryConditions(i);
                
                u.setFValues(obj.updateInitialDisplacement(bc,u));

                residu = 1; residuVec = [];
                obj.elasticity.computeDamageEvolutionParam(u);
               
                Res  = obj.elasticity.computeResidual(u,bc);                
                Dres = obj.elasticity.computeDerivativeResidual(u,bc);
                residu0 = norm(Res);
                
                while (residu >= obj.tolerance)
                    [uNew,uNewVec] = obj.computeU(Dres,Res,u,bc);
                                                          
                    u.setFValues(uNew);
                                      
                    Res  = obj.elasticity.computeResidual(u,bc);                    
                    Dres = obj.elasticity.computeDerivativeResidual(u,bc);
                    % How is Dres dependent on the bc?
                    
                    residu = norm(Res(bc.free_dofs))/residu0;

                    obj.elasticity.computeDamageEvolutionParam(u);

                    residuVec(end+1) = residu/residu0;

                    fprintf('Error: %d | %d \n',residu,uNewVec(6));%.evaluate([0;0]));
                end
                obj.elasticity.setROld(r);
                %fprintf('r  = %d \n',r.evaluate([0;0]));

                data.displacement.value(i)  = obj.boundaryConditions.bcValueSet(i);
                damageDomainFun = obj.elasticity.computeDamage();
                damageFun = damageDomainFun.project('P1D',obj.mesh);
                data.damage.maxValue(i)  = max(damageFun.fValues);
                data.damage.minValue(i)  = min(damageFun.fValues);
                data.reaction(i)  = -obj.computeTotalReaction(Dres,uNewVec);
            end
            data.displacement.field = u;
            data.damage.field = damageFun;
            data.residuVec = residuVec;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = 2;
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material  = cParams.material;
            obj.solverParams = cParams.solver; 
            obj.H = cParams.H;
            obj.r0 = ConstantFunction.create(cParams.r0,obj.mesh);
        end
        
        function bc = updateBoundaryConditions (obj,i)
            bc = obj.boundaryConditions.nextStep(i);
            obj.elasticity.updateBoundaryConditions(bc);
        end

        function createFunctionals(obj)
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.H = obj.H;
            s.r0 = obj.r0;
            
            s.quadOrder = obj.quadOrder;
            s.boundaryConditions = obj.boundaryConditions;

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