classdef ElasticityFunctional < handle
    
    properties (Access = private)
        functionals
        quadOrder
    end
    
    methods (Access = public)
        
        function obj = ElasticityFunctional(cParams)
            obj.init(cParams)
        end

        function Etot = computeCost(obj,u,bc)
            fExt = bc.pointloadFun;
            E    = obj.computeEnergies(u,fExt);
            Etot = sum(E);
        end

        function E = computeEnergies(obj,u,fExt)
            Eint = obj.functionals.intE.computeFunctional(u,obj.quadOrder);
            Wext = obj.functionals.extWork.computeFunctional(u,fExt,obj.quadOrder);
            E = [Eint,Edis,Ereg,Wext];
        end

        function RHS = computeResidual(obj,u,bc)
            fExt = bc.pointloadFun;
            Fint = obj.functionals.intE.computeGradientDisplacement(u,obj.quadOrder);
            Fext = obj.functionals.extWork.computeGradient(u,fExt,obj.quadOrder);
            RHS  = Fint - Fext;
        end

        function LHS = computeHessian(obj,u)
            LHS  = obj.functionals.intE.computeHessianDisplacement(u,obj.quadOrder);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = cParams.quadOrder;
            obj.functionals.extWork = ExternalWorkFunctional(cParams);
            obj.functionals.intE    = NeohookeanFunctional(cParams);
        end

    end

end