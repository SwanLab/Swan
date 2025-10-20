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
                f = x.getFValuesByNode;
                connec = x.mesh.connec;
                dofsC = reshape(connec',1,[]);
                xProj = reshape(f(dofsC,:)',1,[]);
            else
                LHS = obj.computeLHS(x);
                RHS = obj.computeRHS(x);
                xProj = LHS\RHS;
            end
            s.mesh  = x.mesh;
            s.order = obj.order;
            s.ndimf = x.ndimf;
            s.fValues = full(xProj);
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
                    a = repmat(a,1,fun.ndimfTotal);
                    LHS = spdiags(a',0,length(a),length(a));
                otherwise
                    test   = LagrangianFunction.create(fun.mesh, fun.ndimf, obj.order);
                    trial  = test;
                    if length(fun.ndimf) == 1
                        f    = @(u,v) DP(v,u);
                    else
                        f    = @(u,v) DDP(v,u);
                    end
                    LHS = IntegrateLHS(f,test,trial,fun.mesh,'Domain',2);
            end
        end

        function itIs = isP1toP1Dprojection(obj,x)
            itIs = isa(x,'LagrangianFunction') && strcmp(x.order, 'P1') && strcmp(obj.order, 'P1D');
        end

        function RHS = computeRHS(obj,fun)
            test = LagrangianFunction.create(fun.mesh,fun.ndimf,obj.order);
            if length(fun.ndimf) == 1
                f    = @(v) DP(fun,v);
            else
                f    = @(v) DDP(fun,v);
            end
            RHS  = IntegrateRHS(f,test,test.mesh,'Domain',2);
        end

        function ord = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
        end

    end

end