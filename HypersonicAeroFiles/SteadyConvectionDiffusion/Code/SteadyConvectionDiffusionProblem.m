classdef SteadyConvectionDiffusionProblem < handle

    properties (Access = public)
        trial
    end

    properties (Access = private)
        dirValues
        mesh
        stab
    end

    properties (Access = private)
        source
    end

    methods (Access = public)
        function obj = SteadyConvectionDiffusionProblem(cParams)
            obj.init(cParams);
            obj.computeSourceTerm(cParams);
        end

        function compute(obj,a,nu)
            K = obj.computeLHS(a,nu);
            f = obj.computeRHS(a,nu);
            obj.solveSystem(K,f);
        end

        function plot(obj)
            x  = obj.mesh.coord;
            uh = obj.trial.fValues;
            if obj.trial.order == "P1"
                plot(x,uh,'-','LineWidth',1.5)
            else
                [x0,y0]=obj.plotQuadraticElements();
                plot(x0,y0,'-','LineWidth',1.5)
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.dirValues = cParams.dirValues;
            obj.mesh      = cParams.mesh;
            obj.trial     = cParams.trial.copy();
            obj.stab      = cParams.stab;
        end

        function computeSourceTerm(obj,cParams)
            sH         = cParams.sHandle;
            obj.source = AnalyticalFunction.create(sH,1,obj.mesh);
        end

        function LHS = computeLHS(obj,a,nu)            
            s.trial = obj.trial;
            s.mesh  = obj.mesh;
            s.type  = ['ConvDif',obj.stab];
            lhs     = LHSintegrator.create(s);
            LHS     = lhs.compute(a,nu);
        end

        function RHS = computeRHS(obj,a,nu)
            s.mesh = obj.mesh;
            s.type = ['ConvDif',obj.stab];
            int    = RHSintegrator.create(s);
            test   = obj.trial;
            RHS    = int.compute(a,nu,obj.source,test);
        end

        function solveSystem(obj,K,f)
            nDofs   = obj.trial.nDofs;
            bc      = obj.createBoundaryConditions1D();
            dirDofs = bc.dirichlet_dofs;
            b       = bc.dirichlet_vals;
            A       = zeros(2,nDofs);
            A(1,dirDofs(1)) = 1;
            A(2,dirDofs(2)) = 1;
            Ktot = [K A';A zeros(2,2)];
            ftot = [f;b];
            sol  = Ktot\ftot;
            obj.trial.fValues = sol(1:end-2);
        end

        function [x0,y0] = plotQuadraticElements(obj)
            xg     = -1:0.1:1;
            xElg   = obj.mesh.computeXgauss(xg);
            yElg   = obj.trial.evaluate(xg);
            [x0,v] = unique(xElg);
            y0     = yElg(v);
        end

        function bc = createBoundaryConditions1D(obj)
            xMin       = min(obj.mesh.coord(:,1));
            xMax       = max(obj.mesh.coord(:,1));
            isDirLeft  = @(coor)  abs(coor(:,1))==xMin;
            isDirRight = @(coor)  abs(coor(:,1))==xMax;

            sDir{1}.domain    = @(coor) isDirLeft(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = obj.dirValues.left;

            sDir{2}.domain    = @(coor) isDirRight(coor);
            sDir{2}.direction = 1;
            sDir{2}.value     = obj.dirValues.right;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end