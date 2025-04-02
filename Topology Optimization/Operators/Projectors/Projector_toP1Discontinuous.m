classdef Projector_toP1Discontinuous < Projector

    methods (Access = public)

        function obj = Projector_toP1Discontinuous(cParams)
            obj.init(cParams);
        end

       function xP1D = project(obj, x)
            xP1D = P1DiscontinuousFunction.create(obj.mesh,x.ndimf);             
            if isa(x,'LagrangianFunction') && strcmp(x.order, 'P1')
                f = x.fValues;
                connec = obj.mesh.connec;
                dofsC = reshape(connec',1,[]);
                xProj = f(dofsC,:);
            else
                LHS   = obj.computeLHS(xP1D);
                RHS   = obj.computeRHS(xP1D,x);
                xProj = LHS\RHS;
                xProj = reshape(xProj',xP1D.ndimf,[])';
            end
            xP1D.fValues  = xProj;
        end

    end

    methods (Access = private)

        function LHS = computeLHS(obj,test)
            s.type  = 'MassMatrix';
            s.mesh  = test.mesh;
            s.test  = test;
            s.trial = test.copy();
            s.quadratureOrder = 2;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,test,fun)           
            s.quadType = 2;
            s.type = 'ShapeFunction';            
            s.mesh = test.mesh;
            int    = RHSintegrator.create(s);
            RHS    = int.compute(fun,test);
        end

    end

end