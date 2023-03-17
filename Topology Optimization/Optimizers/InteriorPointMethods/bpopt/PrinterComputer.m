classdef PrinterComputer < handle
    properties (Access = public)
              bp
        iter
        x 
        lam
        zL
        zU
        alpha_du
        alpha_pr
        s 
        xL
        xU
        bL
        bU
        screen
        n
        grad
        jacobian
        merit
        objective
        theta
        du
        logmu  
    end
    properties (Access = private)

    end

    methods (Access = public)
        function obj = PrinterComputer(cParams)
            obj.init(cParams);
        end

        function print(obj)
            obj.loadParameters();
            obj.computeObjectiveGradient();
            obj.computeJacobian();
            obj.computedu();
            obj.computeMeritFunction();
            obj.computeObjective();
            obj.computeLogarithmicMu();
            obj.computeTheta();
            obj.printResults();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams.bp;
            obj.iter = cParams.iter;
            obj.x = cParams.x;
            obj.lam = cParams.lam;
            obj.zL = cParams.zL;
            obj.zU = cParams.zU;
            obj.alpha_pr = cParams.alpha_pr;
            obj.alpha_du = cParams.alpha_du;
            obj.s = cParams.s;
            obj.xL = cParams.xL;
            obj.xU = cParams.xU;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
        end

        function loadParameters(obj)
            obj.n = size(obj.x,2);
            obj.screen = 1;
        end

        function computeObjectiveGradient(obj)
            u.x = obj.x;
            u.s = obj.s;
            u.bp = obj.bp;
            gradC = GradientComputer(u);
            gradC.create();
            obj.grad = gradC.objGradient;
        end

        function computeJacobian(obj)
            u.x = obj.x;
            u.s = obj.s;
            u.bp = obj.bp;
            u.bL = obj.bL;
            u.bU = obj.bU;
            jac = JacobianComputer(u);
            jac.compute();
            obj.jacobian = jac.pd;
        end

        function computedu(obj)
            obj.du = sum(abs(obj.grad' + obj.jacobian'*obj.lam' - obj.zL' + obj.zU'));
        end

        function computeMeritFunction(obj)
            u.x = obj.x;
            u.s = obj.s;
            u.xL = obj.xL;
            u.xU = obj.xU;
            u.bL = obj.bL;
            u.bU = obj.bU;
            u.bp = obj.bp;
            me = MeritComputer(u);
            me.compute();
            obj.merit = me.merit;
        end

        function computeObjective(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            object = ObjectiveFunctionComputer(u);
            object.compute();
            obj.objective = object.objectiveFunc;
        end

        function computeLogarithmicMu(obj)
            obj.logmu = log10(obj.bp.mu);
        end

        function computeTheta(obj)
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            u.bp = obj.bp;
            th = ThetaComputer(u);
            th.compute();
            obj.theta = th.theta;
        end

        function printResults(obj)
            if (mod(obj.iter,15)==0)
                fprintf(obj.screen,'   Iter       Merit   Objective   log10(mu)        Pcvg        Dcvg    alpha_pr    alpha_du\n');
            end
            fprintf(obj.screen,'  %5i %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n',obj.iter,obj.merit,obj.objective,obj.logmu,obj.theta,obj.du,obj.alpha_pr,obj.alpha_du);
        end
    end
end