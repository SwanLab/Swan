classdef CohesiveFunctional < handle

    properties (Access = private)
        functionals
        quadOrder
        test
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
                vals = bc.tractionFun.fValues();
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.nnodes,u.mesh.ndim));
            end
            Eint = obj.functionals.energy.computeCost(u);
            Ecoh = obj.functionals.cohesive.computeCost(u,102);
            Wext = obj.functionals.extWork.computeCost(u,fExt,obj.quadOrder);
            Etot = Eint+Ecoh+Wext;
        end

        function RHS = computeGradient(obj,u,bc)
            fExt = bc.tractionFun;
            if ~isempty(bc.tractionFun)
                vals = bc.tractionFun.fValues();
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.nnodes,u.mesh.ndim));
            end
            Fint = obj.functionals.energy.computeGradient(u);
            Fext = obj.functionals.extWork.computeGradient(u,fExt,obj.quadOrder);
            Fcoh = obj.functionals.cohesive.computeResidual(u,102);
            RHS  = Fint-Fext+Fcoh;
        end

        function [KSec,KTan] = computeHessian(obj,u)
            Kelas = obj.functionals.energy.computeHessian();
            [KcohSec, KcohTan]  = obj.functionals.cohesive.computeDerivativeResidual(u,102);
            KTan   = Kelas+KcohTan;
            KSec   = Kelas+KcohSec;
        end

        function updateDamageOld(obj,u)
            obj.functionals.cohesive.updateDamageOld(u);
        end
       
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.quadOrder = cParams.quadOrder;
            obj.test      = cParams.test;
            cParams.testSpace.u = obj.test;
            obj.functionals.extWork  = ExternalWorkFunctional(cParams);
            obj.functionals.energy   = LinearElasticityFunctional(cParams);
            obj.functionals.cohesive = CohesiveTermFunctional(cParams);
        end
        
    end
    
end