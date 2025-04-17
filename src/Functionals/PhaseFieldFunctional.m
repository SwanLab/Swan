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
            E    = obj.computeEnergiesFunctional(u,phi,fExt);
            Etot = sum(E);
        end
        
        function E = computeEnergiesFunctional(obj,u,phi,fExt)
            Eint = obj.functionals.energy.computeFunctional(u,phi,obj.quadOrder);
            Edis = obj.functionals.localDamage.computeFunctional(phi,obj.quadOrder);
            Ereg = obj.functionals.nonLocalDamage.computeFunctional(phi,obj.quadOrder);
            Wext = obj.functionals.extWork.computeFunctional(u,fExt,obj.quadOrder);
            E = [Eint,Edis,Ereg,Wext];
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
                obj.functionals.energy     = PhaseFieldInternalEnergySplitFunctional(cParams);
            else
                obj.functionals.energy     = PhaseFieldInternalEnergyFunctional(cParams);
            end
 

        end
        
    end
    
end