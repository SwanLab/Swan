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
            fExt = bc.tractionFun;
            if ~isempty(bc.tractionFun)
                vals = bc.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.nnodes,u.mesh.ndim));
            end
            E    = obj.computeEnergies(u,fExt);
            Etot = sum(E);
        end

        function E = computeEnergies(obj,u,fExt)
            Eint = obj.functionals.intE.computeCost(u,obj.quadOrder);
            Wext = obj.functionals.extWork.computeCost(u,fExt,obj.quadOrder);
            E = [Eint,Wext];
        end

        function RHS = computeGradient(obj,u,bc)
            fExt = bc.tractionFun;
            if ~isempty(bc.tractionFun)
                vals = bc.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.ndim,u.mesh.nnodes)');
            end
            Fint = obj.functionals.intE.computeGradient(u,obj.quadOrder);
            Fext = obj.functionals.extWork.computeGradient(u,fExt,obj.quadOrder);
            RHS  = Fint - Fext;
        end

        function LHS = computeHessian(obj,u)
            LHS  = obj.functionals.intE.computeHessian(u,obj.quadOrder);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = cParams.quadOrder;
            if strcmp(cParams.matProp.type,'Neohookean')
                cParams.material = cParams.matProp;
                obj.functionals.intE = NeohookeanFunctional(cParams);
            elseif strcmp(cParams.matProp.type,'Elastic')
                cParams.material = cParams.matTensor;
                obj.functionals.intE = LinearElasticityFunctional(cParams);
            end
            obj.functionals.extWork = ExternalWorkFunctional(cParams);
        end

    end

end