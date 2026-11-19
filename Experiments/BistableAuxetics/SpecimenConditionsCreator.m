classdef SpecimenConditionsCreator < handle

    properties (GetAccess = public, SetAccess = private)
        type
        bcValues
        nodeCoordinates
    end

    properties (Access = private)
        mesh
        order
        geometry
        step
        constraints
    end

    methods (Access = public)

        function obj = SpecimenConditionsCreator(mesh,cParams)
            obj.init(mesh,cParams);
            obj.defineConstraints();
        end

        function bc = nextStep(obj)
            obj.step = min(obj.step+1,numel(obj.bcValues));
            bc = obj.computeConditions(obj.bcValues(obj.step));
        end

        function bc = computeConditions(obj,value)
            dirichlet = cell(1,numel(obj.constraints));
            for i = 1:numel(obj.constraints)
                dirichlet{i} = obj.createDirichletCondition(obj.constraints(i),value);
            end
            s.mesh         = obj.mesh;
            s.dirichletFun = [dirichlet{:}];
            s.pointloadFun = [];
            s.periodicFun  = [];
            bc = BoundaryConditions(s);
        end

        function force = computeConjugateForce(obj,fInt)
            force = 0;
            for c = obj.constraints
                force = force + c.factor*sum(fInt(obj.computeDofs(c)));
            end
        end

        function dofs = computeRestrictedDofs(obj)
            dofs = [];
            for c = obj.constraints
                dofs = [dofs;obj.computeDofs(c)];
            end
            dofs = unique(dofs);
        end

    end

    methods (Access = private)

        function init(obj,mesh,cParams)
            obj.mesh     = mesh;
            obj.type     = cParams.type;
            obj.bcValues = cParams.values;
            obj.geometry = cParams.geometry;
            obj.order    = cParams.order;
            obj.step     = 0;
            obj.computeNodeCoordinates();
            obj.constraints = struct('domain',{},'direction',{},'factor',{});
        end

        function defineConstraints(obj)
            g = obj.geometry;
            switch obj.type
                case 'UniaxialRollers'
                    obj.addRollers();
                    obj.lockRigidTranslation();
                case 'ClampedEdges'
                    obj.addClampedEdges();
                case 'RollersFixedSides'
                    obj.addRollers();
                    obj.addFixedSides();
                case 'ClampedFixedSides'
                    obj.addClampedEdges();
                    obj.addFixedSides();
                case 'ClampedWalledSides'
                    obj.addClampedEdges();
                    obj.addSideWalls();
                case 'UniaxialVertices'
                    obj.addVertices(obj.verticesOnGrips(g.bottomGrips),2,@(v) 0);
                    obj.addVertices(obj.verticesOnGrips(g.topGrips),2,@(v) g.height);
                    obj.lockRigidTranslation();
                case 'BiaxialVertices'
                    obj.addAffineVertices(g.latticeVertices);
                case 'ClampedBottomPunch'
                    obj.addClampedBottom();
                    obj.addConstraint(obj.gripDomain(g.topGrips),2,-1);
                case 'BiaxialCorners'
                    v = g.latticeVertices;
                    obj.addAffineVertices(v(convhull(v(:,1),v(:,2)),:));
            end
        end

        function addRollers(obj)
            g = obj.geometry;
            obj.addConstraint(obj.gripDomain(g.bottomGrips),2,0);
            obj.addConstraint(obj.gripDomain(g.topGrips),2,g.height);
        end

        function addClampedEdges(obj)
            g = obj.geometry;
            obj.addRollers();
            obj.addConstraint(obj.gripDomain(g.bottomGrips),1,0);
            obj.addConstraint(obj.gripDomain(g.topGrips),1,0);
        end

        function addClampedBottom(obj)
            bottom = obj.gripDomain(obj.geometry.bottomGrips);
            obj.addConstraint(bottom,1,0);
            obj.addConstraint(bottom,2,0);
        end

        function addFixedSides(obj)
            obj.addConstraint(obj.sideDomain(obj.geometry.sideGrips),1,0);
        end

        function addSideWalls(obj)
            x    = obj.nodeCoordinates(:,1);
            tol  = 1e-8*obj.geometry.sideLength;
            xLim = [min(x),max(x)];
            obj.addConstraint(@(c) abs(c(:,1)-xLim(1)) < tol | abs(c(:,1)-xLim(2)) < tol,1,0);
        end

        function addAffineVertices(obj,vertices)
            c = obj.geometry.center;
            obj.addVertices(vertices,1,@(v) v(1)-c(1));
            obj.addVertices(vertices,2,@(v) v(2)-c(2));
        end

        function addVertices(obj,vertices,direction,factorFun)
            r = 0.4*obj.geometry.cornerLigament;
            for i = 1:size(vertices,1)
                v = vertices(i,:);
                domain = @(x) (x(:,1)-v(1)).^2 + (x(:,2)-v(2)).^2 < r^2;
                obj.addConstraint(domain,direction,factorFun(v));
            end
        end

        function addConstraint(obj,domain,direction,factor)
            assert(any(domain(obj.nodeCoordinates)),'Constraint without nodes');
            c.domain    = domain;
            c.direction = direction;
            c.factor    = factor;
            obj.constraints(end+1) = c;
        end

        function lockRigidTranslation(obj)
            bottom = obj.geometry.bottomGrips;
            target = [mean(bottom(:,1:2),'all'),min(bottom(:,3))];
            x = obj.nodeCoordinates;
            isOnBottom = obj.gripDomain(bottom);
            onBottom = find(isOnBottom(x));
            [~,i] = min(sum((x(onBottom,:)-target).^2,2));
            p = x(onBottom(i),:);
            domain = @(c) abs(c(:,1)-p(1)) + abs(c(:,2)-p(2)) < 1e-10;
            obj.addConstraint(domain,1,0);
        end

        function v = verticesOnGrips(obj,grips)
            v = obj.geometry.latticeVertices;
            isOnGrip = obj.gripDomain(grips);
            v = v(isOnGrip(v),:);
        end

        function domain = gripDomain(obj,grips)
            tol = 1e-8*obj.geometry.sideLength;
            domain = @(x) obj.isOnGrips(x,grips,tol);
        end

        function domain = sideDomain(obj,grips)
            tol = 1e-8*obj.geometry.sideLength;
            domain = @(x) obj.isOnGrips(x(:,[2 1]),grips,tol);
        end

        function d = createDirichletCondition(obj,c,value)
            s.domain    = c.domain;
            s.direction = c.direction;
            s.value     = c.factor*value;
            s.order     = obj.order;
            d = DirichletCondition(obj.mesh,s);
        end

        function computeNodeCoordinates(obj)
            x = LagrangianFunction.create(obj.mesh,2,obj.order).getDofCoord();
            obj.nodeCoordinates = x(1:2:end,:);
        end

        function dofs = computeDofs(obj,c)
            nodes = find(c.domain(obj.nodeCoordinates));
            dofs  = 2*(nodes-1) + c.direction;
        end

    end

    methods (Static, Access = private)

        function is = isOnGrips(x,grips,tol)
            is = false(size(x,1),1);
            for i = 1:size(grips,1)
                inX = x(:,1) > grips(i,1)-tol & x(:,1) < grips(i,2)+tol;
                inY = abs(x(:,2)-grips(i,3)) < tol;
                is  = is | (inX & inY);
            end
        end

    end

end
