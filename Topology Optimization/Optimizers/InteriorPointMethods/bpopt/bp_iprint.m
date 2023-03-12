classdef bp_iprint < handle
    properties (Access = public)
        
    end
    properties (Access = private)
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

    methods (Access = public)
        function obj = bp_iprint(cParams)
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
            gradC = bp_objgrad(obj);
            gradC.compute();
            obj.grad = gradC.objGradient;
        end

        function computeJacobian(obj)
            jac = bp_jac(obj);
            jac.compute();
            obj.jacobian = jac.pd;
        end

        function computedu(obj)
            obj.du = sum(abs(obj.grad' + obj.jacobian'*obj.lam' - obj.zL' + obj.zU'));
        end

        function computeMeritFunction(obj)
            me = bp_merit(obj);
            me.compute();
            obj.merit = me.merit;
        end

        function computeObjective(obj)
            object = bp_obj(obj);
            object.compute();
            obj.objective = object.objectiveFunc;
        end

        function computeLogarithmicMu(obj)
            obj.logmu = log10(obj.bp.mu);
        end

        function computeTheta(obj)
            th = bp_theta(obj);
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