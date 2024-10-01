classdef NonLinearFilterEllipsev2 < handle
    
    properties (Access = private)
        mesh
        trial
        q
        epsilon
    end

    properties (Access = private)
        M1
        M2
        RHSChi
        Kq
        RHS2
    end

    methods (Access = public)
        function obj = NonLinearFilterEllipsev2(cParams)
            obj.init(cParams);
            obj.createMassMatrixFirstDirection();
            obj.createMassMatrixSecondDirection();
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            obj.createRHSChiFirstDirection(fun,quadOrder);
            iter = 1;
            tolerance = 1;
            fr = 0.1;
            while tolerance >= 1e-5 
                oldRho = obj.trial.fValues;
                obj.createKqFirstDirection(quadOrder);
                obj.solveFirstDirection(fr);
                obj.createRHSSecondDirection(quadOrder);
                obj.solveSecondDirection(fr);
                tolerance = norm(obj.trial.fValues - oldRho)/norm(obj.trial.fValues); 
                iter = iter + 1;
                disp(iter);  
                disp(tolerance);
             end
           
           obj.trial.plot
           obj.q.plot
%            obj.differentPlots( chiValues, E1Values, E2Values, rhoValues, ...
%             qValues,iterations);
            xF.fValues  = obj.trial.fValues;
            %disp(iter);
            %disp(tolerance);
            

        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                obj.computeLHS();
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial   = LagrangianFunction.create(cParams.mesh, 1, 'P1'); % rho_eps
            obj.q       = LagrangianFunction.create(cParams.mesh, 2, 'P0'); % 2 = geom dim
            obj.mesh    = cParams.mesh;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end
        function e = updateE1(obj)
                 gRho = Grad(obj.trial);
                 energy_qDRho = gRho.*obj.q;
                 e = Integrator.compute(energy_qDRho,obj.mesh,2);

        end

        function createMassMatrixFirstDirection(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            s.quadratureOrder = 2;
            LHS     = LHSintegrator.create(s);
            obj.M1   = LHS.compute();
        end

        function createMassMatrixSecondDirection(obj)
            s.type  = 'MassMatrix'; % AnisotropicMassMatrix
            s.mesh  = obj.mesh;
            s.test  = obj.q;
            s.trial = obj.q;
            s.quadratureOrder = 2;
            LHS     = LHSintegrator.create(s);
            obj.M2   = LHS.compute();
        end

        function createRHSChiFirstDirection(obj,fun,quadOrder)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(fun,test);
            obj.RHSChi   = rhs;
        end


        function createKqFirstDirection(obj, quadOrder) 
            s.mesh = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            % qRotated = RotatedVector(invA,obj.q);
            rhs        = int.compute(qRotated, test);
            obj.Kq = rhs;
        end

        function createRHSSecondDirection(obj,quadOrder)
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            nablaRho   = Grad(obj.trial);
            %rotatedNablaRho = RotatedVector(invA,nablaRho);
            test       = obj.q;
            rhs        = int.compute(rotatedNablaRho,test);
            obj.RHS2   = -obj.epsilon^2*rhs;
        end



        function solveFirstDirection(obj,fr)
            LHS = obj.M1;
            RHS = obj.RHSChi + obj.Kq;
            rhoi = (1-fr).*obj.trial.fValues + fr.*(LHS\RHS); % hard coded direct solver
            obj.trial.fValues = rhoi;
        end

        function solveSecondDirection(obj,fr)
            LHS = obj.M2;
            RHS = obj.RHS2;
            
            obj.q.fValues = reshape(obj.q.fValues',1,[])';
            qi = (1-fr).*obj.q.fValues +fr.* (LHS \ RHS);
            obj.q.fValues = reshape(qi,[2,obj.mesh.nelem])';
            %obj.q.fValues = reshape(qi,[2,obj.mesh.nnodes])';% 2 = geo dim; P0 q obj.mesh.nelem / P1 q obj.mesh.nnodes
        end


        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end

        function differentPlots(obj,chi, E1,E2,rho,q,iter)
            tiledlayout(3, 2); 
            nexttile;
            plot(iter, E1, 'DisplayName', 'q grad(rho)');
            xlabel('Iteration');
            ylabel('Energy');
            title('Energy Values Over Iterations');
            legend('show'); % Display legend based on 'DisplayName'
            set(gca, 'YScale', 'log')
            grid on;

            % Second subplot (for E2values)
            %subplot(2, 1, 2);  % 2 rows, 1 column, plot 2
            nexttile;
            plot(iter, E2, 'DisplayName', 'q^2');
            xlabel('Iteration');
            ylabel('Energy');
            title('Energy Values Over Iterations');
            legend('show'); % Display legend based on 'DisplayName'
            set(gca, 'YScale', 'log')
            grid on;

             nexttile;
             plot(iter, chi, 'DisplayName', 'Chi');
             xlabel('Iteration');
             ylabel('Energy');
             title('Energy Values Over Iterations');
             legend('show'); % Display legend based on 'DisplayName'
             set(gca, 'YScale', 'log')
             grid on;

            nexttile;
            plot(iter, rho, 'DisplayName', 'rho');
            xlabel('Iteration');
            ylabel('rho');
            title('Rho Values Over Iterations');
            legend('show'); % Display legend based on 'DisplayName'
            set(gca, 'YScale', 'log')
            grid on;

            nexttile;
            plot(iter, q, 'DisplayName', 'q');
            xlabel('Iteration');
            ylabel('q');
            title('q Values Over Iterations');
            legend('show'); % Display legend based on 'DisplayName'
            set(gca, 'YScale', 'log')
            grid on;

        end


    end
end