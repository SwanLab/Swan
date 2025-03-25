classdef PhaseFieldFunctional < handle
    
    properties (Access = private)
        functionals
        quadOrder
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldFunctional(cParams)
            obj.init(cParams)
        end

        function Etot = computeCostFunctional(obj,u,phi,bc)
            fExt = bc.pointloadFun;
            Wext = obj.computeExternalWork(u,fExt);
            Eint = obj.computeInternalEnergy(u,phi);
            Edis = obj.computeDissipationEnergy(phi);
            Ereg = obj.computeRegularisationEnergy(phi);
            Etot = Eint + Edis + Ereg - Wext;
        end
        
        function Wext = computeExternalWork(obj,u,fExt)
            Wext = obj.functionals.extWork.computeFunctional(u,fExt,obj.quadOrder);
        end
        
        function Eint = computeInternalEnergy(obj,u,phi)
            Eint = obj.functionals.energy.computeFunctional(u,phi,obj.quadOrder);
        end
        
        function Edis = computeDissipationEnergy(obj,phi)
            Edis = obj.functionals.localDamage.computeFunctional(phi,obj.quadOrder);
        end
        
        function Ereg = computeRegularisationEnergy(obj,phi)
            Ereg = obj.functionals.nonLocalDamage.computeFunctional(phi,obj.quadOrder);
        end
        
        
        function LHS = computeElasticLHS(obj,u,phi)
            LHS  = obj.functionals.energy.computeHessianDisplacement(u,phi,obj.quadOrder);
        end
        
        function RHS = computeElasticRHS(obj,u,phi,bc)
            fExt = bc.pointloadFun;
            Fint = obj.functionals.energy.computeGradientDisplacement(u,phi,obj.quadOrder);
            Fext = obj.functionals.extWork.computeGradient(u,fExt,obj.quadOrder);
            RHS  = Fint - Fext;
        end
        
        function LHS = computePhaseFieldLHS(obj,u,phi)
            Mi  = obj.functionals.energy.computeHessianDamage(u,phi,obj.quadOrder);
            Md  = obj.functionals.localDamage.computeHessian(phi,obj.quadOrder);
            K   = obj.functionals.nonLocalDamage.computeHessian(phi,obj.quadOrder);
            LHS = Mi + Md + K;
        end
        
        function RHS = computePhaseFieldRHS(obj,u,phi)
            Fi  = obj.functionals.energy.computeGradientDamage(u,phi,obj.quadOrder);
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
            if cParams.energySplit
                obj.functionals.energy     = InternalEnergySplitFunctional(cParams);
            else
                obj.functionals.energy     = InternalEnergyFunctional(cParams);
            end
 

        end
        
    end
    
end