classdef NonLinearFilterEllipsev2 < handle
    
    properties (Access = private)
        mesh
        trial
        q
        epsilon
        A
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
            obj.A = cParams.A;
            obj.createMassMatrixFirstDirection();
            obj.createMassMatrixSecondDirection();  
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            obj.createRHSChiFirstDirection(fun,quadOrder);
            iter = 1;
            tolerance = 1;
            fr = 0.01;
%             filename = 'Ellipse90.gif'; 
            while tolerance >= 1e-5
                oldRho = obj.trial.fValues;
                obj.createKqFirstDirection(quadOrder);
                obj.solveFirstDirection(fr);
                obj.createRHSSecondDirection(quadOrder);
                obj.solveSecondDirection(fr);
                tolerance = norm(obj.trial.fValues - oldRho)/norm(obj.trial.fValues); 
                disp(tolerance);
%                 if mod(iter, 5) == 0
%                     % Create the plot without displaying it
%                     figure('Visible', 'off');  % Create a figure without displaying it
%                     obj.trial.plot;            % Call the plotting method
%                     title(['Iteration: ', num2str(iter)]);
%                     drawnow;                   % Update the plot immediately
%                     % Capture the plot as an image      
%                     frame = getframe(gcf);       % Capture the current figure frame
%                     im = frame2im(frame);        % Convert the frame to an image
%                     [imind, cm] = rgb2ind(im, 256); % Convert the image to an indexed image
%             
%                     if iter == 5  % Check if it's the first frame to create GIF
%                         imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
%                     else
%                         % Append subsequent frames to the GIF file
%                         imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
%                     end
%             
%                     % Close the figure after capturing
%                     close(gcf);  % Close the figure to avoid displaying it
%                 end
        
        % Increment the iteration counter
                iter = iter + 1;
                disp(iter);
            end
            
            obj.trial.plot
            xF.fValues  = obj.trial.fValues;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                %obj.computeLHS();
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial   = LagrangianFunction.create(cParams.mesh, 1, 'P1'); % rho_eps
            obj.q       = LagrangianFunction.create(cParams.mesh, 2, 'P0'); % 2 = geom dim
            obj.mesh    = cParams.mesh;
            obj.epsilon = 5*cParams.mesh.computeMeanCellSize();
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
            s.type  = 'AnisotropicMassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.q;
            s.trial = obj.q;
            s.quadratureOrder = 2;
            s.A     = obj.A; % A from cParams
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
            qRotated = RotatedVector(obj.A,obj.q);
            rhs        = int.compute(qRotated, test);
            obj.Kq = rhs;
        end

        function createRHSSecondDirection(obj,quadOrder)
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            nablaRho   = Grad(obj.trial);
            rotatedNablaRho = RotatedVector(obj.A,nablaRho);
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
    end
end