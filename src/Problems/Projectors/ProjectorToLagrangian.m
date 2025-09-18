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
            s.fValues = full(xProj);
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
                    test   = LagrangianFunction.create(fun.mesh, fun.ndimf, obj.order);
                    trial  = LagrangianFunction.create(fun.mesh, fun.ndimf, obj.order);
                    f = @(u,v) DP(v,u);
                    LHS = IntegrateLHS(f,test,trial,fun.mesh,2);
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
                    s.quadType = ord;
                    int        = RHSIntegrator.create(s);
                    test       = LagrangianFunction.create(fun.mesh,fun.ndimf,obj.order);
                    RHS        = int.compute(fun,test);
                otherwise
                    test       = LagrangianFunction.create(fun.mesh,fun.ndimf,obj.order);
                    f = @(v) DP(v,fun);
                    RHS = IntegrateRHS(f,test,test.mesh,2);   
            end

        end

        function ord = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
        end

    end

end