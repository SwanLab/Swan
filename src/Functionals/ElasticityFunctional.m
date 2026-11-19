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
            Wext = 0;
            if ~isempty(fExt)
                Wext = obj.functionals.extWork.computeCost(u,fExt,obj.quadOrder);
            end
            E = [Eint,Wext];
            if obj.hasContact()
                E(end+1) = obj.functionals.contact.computeCost(u);
            end
        end

        function RHS = computeGradient(obj,u,bc)
            fExt = bc.tractionFun;
            if ~isempty(bc.tractionFun)
                vals = bc.tractionFun.computeRHS([]);
                fExt = LagrangianFunction.create(u.mesh, u.mesh.ndim,'P1');
                fExt.setFValues(reshape(vals,u.mesh.ndim,u.mesh.nnodes)');
            end
            Fint = obj.functionals.intE.computeGradient(u,obj.quadOrder);
            RHS  = Fint;
            if ~isempty(fExt)
                RHS = RHS - obj.functionals.extWork.computeGradient(u,fExt,obj.quadOrder);
            end
            if obj.hasContact()
                RHS = RHS + obj.functionals.contact.computeGradient(u);
            end
        end

        function resetParameters(obj)
            for f = obj.functionalsWithParameters()
                f{1}.resetParameters();
            end
        end

        function itHas = updateParameters(obj,u)
            itHas = true;
            for f = obj.functionalsWithParameters()
                itHas = f{1}.updateParameters(u) && itHas;
            end
        end

        function LHS = computeHessian(obj,u)
            LHS  = obj.functionals.intE.computeHessian(u,obj.quadOrder);
            if obj.hasContact()
                LHS = LHS + obj.functionals.contact.computeHessian(u);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = cParams.quadOrder;
            if strcmp(cParams.matProp.type,'Neohookean')
                cParams.material = cParams.matProp;
                obj.functionals.intE = NeohookeanFunctional(cParams);
            elseif strcmp(cParams.matProp.type,'ThirdMedium')
                cParams.material = cParams.matProp;
                obj.functionals.intE = ThirdMediumFunctional(cParams);
            elseif strcmp(cParams.matProp.type,'Elastic')
                cParams.material = cParams.matTensor;
                obj.functionals.intE = LinearElasticityFunctional(cParams);
            end
            if strcmp(cParams.testSpace.u.order,'P1')
                obj.functionals.extWork = ExternalWorkFunctional(cParams);
            end
            if isfield(cParams,'contact')
                s = cParams.contact;
                s.mesh = cParams.mesh;
                switch s.type
                    case 'penalty'
                        obj.functionals.contact = SelfContactFunctional(s);
                    case 'barrier'
                        obj.functionals.contact = BarrierContactFunctional(s);
                    case 'augmentedLagrangian'
                        obj.functionals.contact = AugmentedLagrangianContactFunctional(s);
                end
            end
        end

        function f = functionalsWithParameters(obj)
            f = struct2cell(obj.functionals)';
            f = f(cellfun(@(g) ismethod(g,'updateParameters'),f));
        end

        function itHas = hasContact(obj)
            itHas = isfield(obj.functionals,'contact');
        end

    end

end