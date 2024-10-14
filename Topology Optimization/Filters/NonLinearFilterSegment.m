classdef NonLinearFilterSegment < handle
    
    properties (Access = private)
        mesh
        trial
        sVar
        k
        alpha
        beta
    end

    properties (Access = private)
        M
        IntChi
    end

    methods (Access = public)
        function obj = NonLinearFilterSegment(cParams)
            obj.init(cParams);
            obj.createMassMatrix();
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            obj.createRHSChi(fun,quadOrder);
            iter = 1;
            tolerance = 1;
            %fr = 0.1;
            while tolerance >= 1e-5 
                oldRho = obj.trial.fValues;
                % Compute remaining LHSs and RHSs
                % Solve for rhoEps
                % Update s
                tolerance = norm(obj.trial.fValues - oldRho)/norm(obj.trial.fValues); 
                iter = iter + 1;
                disp(iter);  
                disp(tolerance);
             end
           
           obj.trial.plot
            

        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial = LagrangianFunction.create(cParams.mesh, 1, 'P1'); % rho_eps
            obj.sVar  = 0;
            obj.mesh  = cParams.mesh;
            obj.k     = cParams.k;
            obj.alpha = cParams.alpha;
            obj.beta  = cParams.beta;
        end

        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            %s.quadratureOrder = 2;
            LHS     = LHSintegrator.create(s);
            obj.M   = LHS.compute();
        end

        function createRHSChi(obj,fun,quadOrder)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(fun,test);
            obj.IntChi = rhs;
        end
    end
end