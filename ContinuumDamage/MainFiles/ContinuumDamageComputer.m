classdef ContinuumDamageComputer < handle

    properties (Access = private)
        tolerance

        mesh
        boundaryConditions
        material
        %H
        solverParams
        quadOrder 

 
    end

    properties (Access = private)
        tau = 5e-4;
        r0
        r1

        elasticFun
        externalWorkFun
        elasticity
        qLaw

        limIter = 100;
    end

    methods (Access = public)
        function obj = ContinuumDamageComputer(cParams)
            obj.init(cParams)
            obj.defineRfunction();
            obj.defineFunctional()
        end

        function data = compute(obj)
            uFun = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.elasticity.setTestFunctions(uFun);
            data = {};

            nSteps = length(obj.boundaryConditions.bcValueSet);
            for i = 1:nSteps
                fprintf('Step: %d / %d \n',i,nSteps);
                bc = obj.updateBoundaryConditions(i);
                uFun.setFValues(obj.updateInitialDisplacement(bc,uFun));
                isLoading = obj.loadState(i);

                err = 1; iter = 0;
                while (err >= obj.tolerance && iter < obj.limIter)    
                    [res,K,uVec,uFun] = obj.solveU(uFun,bc,isLoading);  

                    err = norm(res);
                    fprintf('Error: %d \n',err);

                    iter = iter+1;
                end
                if (iter >= obj.limIter)
                    fprintf (2,'NOT CONVERGED FOR STEP %d\n',i);
                end
                obj.elasticity.setROld(); 
                [data,dmgFun,rFun,qFun] = obj.getData(data,i,K,uVec,uFun,bc);
                
            end
            data.displacement.field = uFun;
            data.damage.field = dmgFun;
            data.r.field = rFun;
            data.q.field = qFun;

        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.material  = cParams.material;
            obj.solverParams = cParams.solver; 
            % obj.H = cParams.H;
            obj.tolerance = cParams.tol;
            obj.quadOrder = 2;
            obj.qLaw = cParams.damageLaw;
        end

        function defineRfunction(obj)
            obj.r0 = LagrangianFunction.create(obj.mesh,1,'P0');
            fV = obj.qLaw.r0*ones(size(obj.r0.fValues));
            obj.r0.setFValues(fV);

            obj.r1 = LagrangianFunction.create(obj.mesh,1,'P0');
            fV = obj.qLaw.r1*ones(size(obj.r1.fValues));
            obj.r1.setFValues(fV);
        end
       
        function defineFunctional(obj)
            s.mesh     = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.material = obj.material;
            % s.H = obj.qLaw.H;
            s.r0 = obj.r0;
            s.r1 = obj.r1;
            s.quadOrder = obj.quadOrder;
            s.qLaw = obj.qLaw;
            obj.elasticity = ShFunc_ContinuumDamage(s);
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

        function [uOut,uOutVec] = computeDisplacement(obj,LHS,RHS,uIn,bc)            
            uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
            uOutVec = uInVec;

            uInFree = uInVec(bc.free_dofs);
            uOutFree = obj.updateWithNewton(LHS,RHS,uInFree);
           % uOutFree = obj.updateWithGradient(RHS,uInFree);
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

        function isLoading = loadState(obj,i)
           isLoading = i <= obj.boundaryConditions.LoadingBcLength ;
        end
        
        function [data,dmgFun,rFun,qFun] = getData(obj,data,i,K,uVec,uFun,bc)
            data.displacement.value(i)  = obj.boundaryConditions.bcValueSet(i);
            dmgDomainFun = obj.elasticity.getDamage();
            dmgFun = dmgDomainFun.project('P1D');
            data.damage.maxValue(i)  = max(dmgFun.fValues);
            data.damage.minValue(i)  = min(dmgFun.fValues);
           
            rDomainFun = obj.elasticity.getR();
            rFun = rDomainFun.project('P0');
            data.r.maxValue(i) = max(rFun.fValues);
            data.r.minValue(i) = min(rFun.fValues);
            
            qDomainFun = obj.elasticity.getQ();
            qFun = qDomainFun.project('P1D');
            data.q.maxValue(i) = max(qFun.fValues);
            data.q.minValue(i) = min(qFun.fValues);
            
            data.reaction(i)  = obj.computeTotalReaction(K,uVec);
            [data.totalEnergy(i),data.damagedMaterial(i)] = obj.elasticity.computeEnergy(uFun,bc);
            
        end
        function [res,K,uVec,uFun] = solveU(obj,uFun,bc,isLoading)  
            obj.elasticity.computeDamageEvolutionParam(uFun);         
            [res]  = obj.elasticity.computeResidual(uFun,bc);
            [K,resDeriv] = obj.elasticity.computeDerivativeResidual(uFun,bc,isLoading);
            [ufV,uVec] = obj.computeDisplacement(resDeriv,res,uFun,bc);
            uFun.setFValues(ufV);
            [res]  = obj.elasticity.computeResidual(uFun,bc);
        end

        function totReact = computeTotalReaction(obj,LHS,u)
            F = LHS*u;
            isInDown =  (abs(obj.mesh.coord(:,2) - min(obj.mesh.coord(:,2)))< 1e-12);
            nodes = 1:obj.mesh.nnodes;
            totReact = -sum(F(2*nodes(isInDown)));
        end
    end
end