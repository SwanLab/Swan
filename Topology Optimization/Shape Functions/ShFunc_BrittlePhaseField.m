classdef ShFunc_BrittlePhaseField < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        functionals
        quadOrder
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ShFunc_BrittlePhaseField(cParams)
            obj.init(cParams)
        end
        
        %%% Energies %%%

        function Etot = computeCostFunction(obj,u,phi,bc)
            fExt = bc.pointloadFun;
            Wext = obj.computeExternalWork(u,fExt);
            Eint = obj.computeInternalEnergy(u,phi);
            Edis = obj.computeDissipationEnergy(phi);
            Ereg = obj.computeRegularisationEnergy(phi);
            Etot = Eint+Edis+Ereg-Wext;
        end
        
        function Wext = computeExternalWork(obj,u,fExt)
            Wext = obj.functionals.extWork.computeFunction(u,fExt,obj.quadOrder);
        end
        
        function Eint = computeInternalEnergy(obj,u,phi)
            Eint = obj.functionals.energy.computeFunction(u,phi,obj.quadOrder);
        end
        
        function Edis = computeDissipationEnergy(obj,phi)
            Edis = obj.functionals.localDamage.computeFunction(phi,obj.quadOrder);
        end
        
        function Ereg = computeRegularisationEnergy(obj,phi)
            Ereg = obj.functionals.nonLocalDamage.computeFunction(phi,obj.quadOrder);
        end
        
        %%% Derivatives %%%
        
        function LHS = computeElasticLHS(obj,u,phi)
            LHS = obj.functionals.energy.computeHessianDisplacement(u,phi,obj.quadOrder);
        end
        
        function RHS = computeElasticRHS(obj,u,phi,bc)
            fExt     = bc.pointloadFun;
            Fint     = obj.functionals.energy.computeGradientDisplacement(u,phi,obj.quadOrder);
            Fext     = obj.functionals.extWork.computeGradient(u,fExt,obj.quadOrder);
            RHS      = Fint + Fext;
        end
        
        function LHS = computePhaseFieldLHS(obj,u,phi)
            Mi = obj.functionals.energy.computeHessianDamage(u,phi,obj.quadOrder);
            Md = obj.functionals.localDamage.computeHessian(phi,obj.quadOrder);
            K  = obj.functionals.nonLocalDamage.computeHessian(phi,obj.quadOrder);
            LHS = Mi + Md + K;
        end
        
        function res = computePhaseFieldRHS(obj,u,phi)
            Fi = obj.functionals.energy.computeGradientDamage(u,phi,obj.quadOrder);
            Fd = obj.functionals.localDamage.computeGradient(phi,obj.quadOrder); 
            DF = obj.functionals.nonLocalDamage.computeGradient(phi,obj.quadOrder);        
            res = Fi + Fd + DF;
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.quadOrder = cParams.quadOrder;
            obj.functionals.energy         = ShFunc_InternalEnergy(cParams);
            %obj.functional.energy2         = ShFunc_InternalEnergySplit(cParams);
            obj.functionals.localDamage    = ShFunc_LocalDamage(cParams);
            obj.functionals.nonLocalDamage = ShFunc_NonLocalDamage(cParams);
            obj.functionals.extWork        = ShFunc_ExternalWork(cParams);
        end
        
    end
    
end