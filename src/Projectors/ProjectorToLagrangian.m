classdef ProjectorToLagrangian < Projector
    
    properties (Access = private)
        order
    end
    
    methods (Access = public)

        function obj = ProjectorToLagrangian(cParams)
            obj.order = cParams.projectorType;
        end

        function xFun = project(obj, x)
            
            if obj.isP1toP1Dprojection(x)
                f = x.fValues;
                connec = x.mesh.connec;
                dofsC = reshape(connec',1,[]);
                xProj = f(dofsC,:);
            else
                LHS = obj.computeLHS(x);
                RHS = obj.computeRHS(x);
                xProj = LHS\RHS;
                xProj = reshape(xProj,[x.ndimf,numel(xProj)/x.ndimf])';
            end
            s.mesh    = x.mesh;
            s.fValues = xProj;
            s.order = obj.order;
            xFun = LagrangianFunction(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj,fun)
            switch obj.order
                case 'P0'
                    quad = Quadrature.create(fun.mesh,1);
                    dv = fun.mesh.computeDvolume(quad);
                    a = sum(dv(1,:),1);
                    a = repmat(a,1,fun.ndimf);
                    LHS = spdiags(a',0,length(a),length(a));
                otherwise
                    s.mesh  = fun.mesh;
                    s.test  = LagrangianFunction.create(fun.mesh, fun.ndimf, obj.order);
                    s.trial = LagrangianFunction.create(fun.mesh, fun.ndimf, obj.order);
                    s.quadratureOrder = 2; % no
                    s.type  = 'MassMatrix';
                    lhs = LHSintegrator.create(s);
                    LHS = lhs.compute();
            end
        end

        function itIs = isP1toP1Dprojection(obj,x)
            itIs = isa(x,'LagrangianFunction') && strcmp(x.order, 'P1') && strcmp(obj.order, 'P1D');
        end

        function RHS = computeRHS(obj,fun)
            ord = obj.createRHSQuadrature(fun);
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = fun.unfittedMesh;
                    s.type = 'Unfitted';
                otherwise
                    s.mesh = fun.mesh;
                    s.type = 'ShapeFunction';
            end
            s.quadType = ord;
            int        = RHSintegrator.create(s);
            test       = LagrangianFunction.create(fun.mesh,fun.ndimf,obj.order);
            RHS        = int.compute(fun,test);
        end

        function ord = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
        end

    end

end