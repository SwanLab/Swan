
classdef CohesiveFunctional < handle
    

    properties (Access = private)
        functionals
        quadOrder
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveFunctional(cParams)
            obj.init(cParams)
        end

        function Etot = computeCost(obj,u,bc)
            fExt = bc.tractionFun;
            if ~isempty(bc.tractionFun)
                vals = bc.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.nnodes,u.mesh.ndim));
            end
            Eint = obj.functionals.energy.computeCost(u,obj.quadOrder);
            Ecoh = obj.functinonals.cohesive.computeCost(u,obj.quadOrder, ????);
            Wext = obj.functionals.extWork.computeCost(u,fExt,obj.quadOrder);
            Etot = Eint+Ecoh+Wext;
        end

        function RHS = computeGradient(obj,u,bc)
            fExt = bc.tractionFun;
            if ~isempty(bc.tractionFun)
                vals = bc.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.nnodes,u.mesh.ndim));
            end
            Fint = obj.functionals.energy.computeGradient(u,obj.quadOrder);
            Fext = obj.functionals.extWork.computeGradient(u,fExt,obj.quadOrder);
            Fcoh = obj.functionals.cohesive.computeResidual(u,obj.quadOrder);
            RHS  = Fint-Fext+Fcoh;
        end

        function LHS = computeHessian(obj,u)
            Kelas = obj.functionals.energy.computeHessian(u,obj.quadOrder);
            Kcoh  = obj.functionals.cohesive.computeDerivativeResidual(u,obj.quadOrder);
            LHS   = Kelas+Kcoh;
        end
       
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.quadOrder = cParams.quadOrder;
            obj.functionals.extWork  = ExternalWorkFunctional(cParams);
            obj.functionals.energy   = LinearElasticityFunctional(cParams);
            obj.functionals.cohesive = CohesiveTermFunctional(cParams);
        end
        
    end
    
end