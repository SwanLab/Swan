classdef PhaseFieldFunctional < handle
    
    properties (Access = private)
        functionals
        quadOrder
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldFunctional(cParams)
            obj.init(cParams)
        end

        function Etot = computeCost(obj,u,theta,phi,bc)
            fExt = bc.tractionFun;
            if ~isempty(bc.tractionFun)
                vals = bc.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.ndim,u.mesh.nnodes)');
            end
            E    = obj.computeEnergies(u,theta,phi,fExt);
            Etot = sum(E);
        end
        
        function E = computeEnergies(obj,u,theta,phi,fExt)
            Eint = obj.functionals.energy.computeCost(u,theta,phi,obj.quadOrder);
            Edis = obj.functionals.localDamage.computeCost(phi,obj.quadOrder);
            Ereg = obj.functionals.nonLocalDamage.computeCost(phi,obj.quadOrder);
            Wext = obj.functionals.extWork.computeCost(u,fExt,obj.quadOrder);
            E = [Eint,Edis,Ereg,Wext];
        end
        
        
        function LHS = computeElasticLHS(obj,u,theta,phi)
            LHS  = obj.functionals.energy.computeHessianDisplacement(u,theta,phi,obj.quadOrder);
        end
        
        function RHS = computeElasticRHS(obj,u,theta,phi,bc)
            fExt = bc.tractionFun;
            if ~isempty(bc.tractionFun)
                vals = bc.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.ndim,u.mesh.nnodes)');
            end
            Fint = obj.functionals.energy.computeGradientDisplacement(u,theta,phi,obj.quadOrder);
            Fext = obj.functionals.extWork.computeGradient(u,fExt,obj.quadOrder);
            RHS  = Fint - Fext;
        end
        
        function LHS = computePhaseFieldLHS(obj,u,theta,phi)
            Mi  = obj.functionals.energy.computeHessianDamage(u,theta,phi,obj.quadOrder);
            Md  = obj.functionals.localDamage.computeHessian(phi,obj.quadOrder);
            K   = obj.functionals.nonLocalDamage.computeHessian(phi,obj.quadOrder);
            LHS = Mi + Md + K;
        end
        
        function RHS = computePhaseFieldRHS(obj,u,theta,phi)
            Fi  = obj.functionals.energy.computeGradientDamage(u,theta,phi,obj.quadOrder);
            Fd  = obj.functionals.localDamage.computeGradient(phi,obj.quadOrder); 
            DF  = obj.functionals.nonLocalDamage.computeGradient(phi,obj.quadOrder);        
            RHS = Fi + Fd + DF;
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.quadOrder = cParams.quadOrder;
            obj.functionals.localDamage    = LocalDamageFunctional(cParams);
            obj.functionals.nonLocalDamage = NonLocalDamageFunctional(cParams);
            obj.functionals.extWork        = ExternalWorkFunctional(cParams);
            obj.functionals.energy         = PhaseFieldInternalEnergyRotationFunctional(cParams);
        end
        
    end
    
end