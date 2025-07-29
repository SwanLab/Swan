classdef HyperelasticityComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        functional
    end

    properties (Access = private)
        monitor
        updater
    end

    methods (Access = public)

        function obj = HyperelasticityComputer(cParams)
            obj.init(cParams);
            obj.setMonitoring(cParams);
            obj.setOptimizer(cParams)
        end

        function outputData = compute(obj)
            u = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            cost = 0;

            nSteps = length(obj.boundaryConditions.bcValues);
            for iStep = 1:nSteps
                obj.monitor.printStep(iStep,nSteps)
                [u,bc] = obj.updateBoundaryConditions(u);
                [u,F,cost,iterMax] = obj.updater.update(u,bc,cost);
                %[Evec,totE,totF,uBC] = obj.postprocess(iStep,u,phi,F,bc);
                %obj.printAndSave(iStep,totF,uBC,u,phi,Evec,totE,iterMax,cost,tauArray);
            end
            outputData = obj.monitor.data;

            % u  = obj.uFun;
            % iter = 1; nSteps = length(obj.boundaryConditions.bcValues);
            % for iStep = 1:nSteps
            %     [u,bc] = obj.updateBoundaryConditions(u);
            % 
            %     %intEnergy(iter) = obj.neohookeanFun.compute(u);
            %     intEnergy(iter) = obj.linearElasticityFun.compute(u);
            %     Res = obj.computeResidual(u);
            %     resi(iter) = norm(Res);
            % 
            %     hasNotConverged = 1;
            %     while hasNotConverged
            %         uVal = obj.reshapeToVector(u.fValues);
            %         [hess, KR] = obj.computeHessian(u);
            %         Res  = obj.computeResidual(u);
            %         incU = hess\(-Res);
            %         uVal(obj.freeDofs) = uVal(obj.freeDofs) + incU;
            %         u.setFValues(obj.reshapeToMatrix(uVal));
            %         reacFun = obj.computeReactions(KR,u);
            % 
            %         iter = iter+1;
            %         intEnergy(iter) = obj.neohookeanFun.compute(u);
            %         Res        = obj.computeResidual(u);
            %         resi(iter) = norm(Res);
            % 
            %         hasNotConverged = obj.hasNotConverged(intEnergy);
            %     end
            %     normDisp(iStep) = Norm(u,'L2');
            %     normReac(iStep) = Norm(reacFun,'L2');
            %     obj.printFile(iStep, u, reacFun);
            %     nIterPerStep(iStep) = obj.computeNumberOfIterations(iter,nIterPerStep,iStep);
            %     obj.plotStep(intEnergy,nIterPerStep,iStep, normDisp,normReac);
            %     obj.rFun = reacFun;
            % end
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.functional         = cParams.functional;
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay = cParams.monitoring.set;
            s.shallPrint   = cParams.monitoring.print;
            obj.monitor = HyperelasticityMonitoring(s);
        end

        function setOptimizer(obj,cParams)
            s.functional  = obj.functional;
            s.monitor     = obj.monitor;
            s.tolerance   = cParams.tolerance;
            s.maxIter     = cParams.maxIter;
            obj.updater = DisplacementUpdater(s);
        end

        function [u,bc] = updateBoundaryConditions(obj,u)
            bc = obj.boundaryConditions.nextStep();
            u.setFValues(obj.updateInitialDisplacement(u,bc));
        end

        function u = updateInitialDisplacement(~,uOld,bc)
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







        function printFile(obj,iStep, u, reac)
            if obj.print
                fun = {u, reac};
                funNames = {'Displacement', 'Reactions'};
                a.mesh     = obj.mesh;
                a.filename = ['SIM_',obj.fileNameOut,'_',int2str(iStep)];
                a.fun      = fun;
                a.funNames = funNames;
                a.type     = 'Paraview';
                pst = FunctionPrinter.create(a);
                pst.print();
            end
        end

        function reacFun = computeReactions(obj,KR,u)
            uVal = obj.reshapeToVector(u.fValues);
            reactions = zeros(size(uVal));
            reactions(obj.constrainedDofs,1) = -KR*uVal;
            s.fValues = obj.reshapeToMatrix(reactions);
            s.order = 'P1';
            s.mesh  = obj.mesh;
            reacFun = LagrangianFunction(s);
        end
        % 
        % function plotStep(obj,intEnergy,nIterPerStep,iStep, normDisp, normReac)
        % 
        %     f = figure(1);
        %     clf(f)
        %     subplot(1,3,1)
        %     plot(intEnergy)
        %     title('Energies')
        %     xlabel('iteration')
        %     ylabel('Energy')
        %     subplot(1,3,2)
        %     plot(normDisp, normReac)
        %     title('Displacement-reaction')
        %     xlabel('displacement')
        %     ylabel('reactions')
        %     hold on
        %     subplot(1,3,3)
        %     bar(1:iStep, nIterPerStep)
        %     title('Number of iterations to converge $\Delta$ F')
        %     xlabel('step')
        %     ylabel('num. iterations per step')
        %     hold on
        %     drawnow
        % end

        function hasNot = hasNotConverged(obj,energy)
            TOL = 10e-12;
            energyOld = energy(end-1);
            energy    = energy(end);
            hasNot = abs(energy-energyOld)/energy > TOL;
        end



        function u = computeInitialElasticGuess(obj, bc)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material;
            s.dim = '2D';
            s.boundaryConditions = bc;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            fem.solve();
            %             u = obj.reshapeToVector(fem.uFun.fValues);
            u = fem.uFun;
        end





        function [Hff,Hr] = computeHessian(obj,u)
            H = obj.linearElasticityFun.computeHessian(u);
            %H = obj.neohookeanFun.computeHessian(u);
            f = obj.freeDofs;
            r = obj.constrainedDofs;
            Hff = H(f,f);
            Hr  = H(r,:);
        end

        function Rf = computeResidual(obj,u)
            Fext          = obj.computeExternalForces();
            [Fint,FintEl] = obj.computeInternalForces(u);
            R = Fint - Fext;% - react;
            f = obj.freeDofs;
            Rf = R(f);
        end

        function nIter = computeNumberOfIterations(obj,iter,iterOld,iStep)
            if iStep == 1
                nIter = iter-1;
            else
                nIter = iter-1-sum(iterOld(1:(iStep-1)));
            end
        end

        function computeFreeDofs(obj)
            bc = obj.createBoundaryConditions(1);
            nDofs = obj.uFun.nDofs;
            obj.freeDofs = setdiff(1:nDofs,bc.dirichlet_dofs);
            obj.constrainedDofs = bc.dirichlet_dofs;
        end



        function u = applyDirichletToUFun(obj,u,bc)
            u_k = reshape(u.fValues',[u.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;
            u.setFValues(reshape(u_k,[obj.mesh.ndim,obj.mesh.nnodes])');
        end

        function [Fext] = computeExternalForces(obj,perc)
            %             Fext = perc*obj.FextInitial;
            %             Fext = obj.reshapeToVector(Fext);
            Fext = zeros(obj.uFun.nDofs,1);
        end

        function [intFor,intForel] = computeInternalForces(obj,u)
            %intFor = obj.neohookeanFun.computeGradient(u);
            intFor = obj.linearElasticityFun.computeGradient(u);
        end

        function rshp = reshapeToVector(obj, A)
            rshp = reshape(A',[obj.uFun.nDofs,1]);
        end

        function rshp = reshapeToMatrix(obj, A)
            rshp = reshape(A,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

        function plotVector(obj, vect)
            s.fValues = obj.reshapeToMatrix(vect);
            s.mesh = obj.mesh;
            s.order = 'P1';
            a = LagrangianFunction(s);
            a.plot;
        end


        % Other useful stuff later on
        function lambdas = computeStretches(obj)
            GradU2 = ActualGrad(obj.uFun);
            Id = Identity(obj.uFun);
            F = GradU2 + Id;
            C = F'*F;
            lambdas = (Eigen(C).^0.5);
        end

        function sigma = computeCauchyStress(obj)
            mu = obj.material.mu;
            lambda = obj.material.lambda;
            GradU2 = ActualGrad(obj.uFun);
            Id = Identity(obj.uFun);
            F = GradU2 + Id;
            b = F*F';
            jac = Det(F);
            sigma = mu.*(b-Id)./jac + lambda.*(log(jac)).*Id./jac;
            sigma.ndimf = [2 2];
        end


    end

end
