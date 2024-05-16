classdef StokesFlowCylinder < handle

    properties (Access = private)
        backgroundMesh
        center
        radius
        velocityOrder
        pressureOrder
    end

    properties (Access = private)
        innerMesh
        vFun
        pFun
        filter
        material
    end

    methods (Access = public)
        function obj = StokesFlowCylinder(cParams)
            obj.init(cParams);
            obj.createInnerMesh();
            obj.createTrialFunctions();
            obj.createPressureSmoother();
            obj.createMaterial();
        end

        function [v,p] = solve(obj,vInlet,vWalls,dtime)
            dirVel            = obj.computeDirichletConditionsHandles(vInlet,vWalls);
            [dir,dirDofs]     = obj.computeDirichletBCMatrices(dirVel);
            LHS               = obj.computeLHS(dtime);
            RHS               = obj.computeRHS(LHS,dir,dirDofs);
            [vValues,pValues] = obj.solveSystem(LHS,RHS,dir,dirDofs);
            [v,p]             = obj.fillTrialFunctions(vValues,pValues);
        end

        function [pCyl,bMesh] = obtainPressureDistributionAtCylinder(obj,pressureFun)
            [~,~,~,isCyl] = obj.computeBoundariesHandles();
            mesh          = obj.innerMesh;
            nodesCyl      = pressureFun.getDofsFromCondition(isCyl);
            pCylVals      = pressureFun.fValues(nodesCyl,1);
            mesh.computeEdges();
            e            = mesh.edges.nodesInEdges;
            bE = ismember(e,nodesCyl);
            bE = find(prod(bE,2));
            connec = e(bE,:);
            s.coord      = mesh.coord;
            s.connec     = connec;
            s.kFace      = -1;
            bMesh        = Mesh.create(s);
            bMesh        = bMesh.computeCanonicalMesh();
            pCyl         = LagrangianFunction.create(bMesh,1,pressureFun.order);
            pCyl.fValues = pCylVals;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.center         = cParams.center;
            obj.radius         = cParams.radius;
            obj.velocityOrder  = cParams.velocityOrder;
            obj.pressureOrder  = cParams.pressureOrder;
        end

        function createInnerMesh(obj)
            m         = obj.backgroundMesh;
            xpos      = obj.center(1);
            ypos      = obj.center(2);
            r         = obj.radius;
            s.type    = 'Given';
            s.fHandle = @(x) -((x(1,:,:)-xpos).^2+(x(2,:,:)-ypos).^2-r^2);
            g         = GeometricalFunction(s);
            lsFun     = g.computeLevelSetFunction(m);

            sUm.backgroundMesh = m;
            sUm.boundaryMesh   = m.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(lsFun.fValues);
            obj.innerMesh = uMesh.createInnerMesh();
        end

        function createTrialFunctions(obj)
            mesh     = obj.innerMesh;
            obj.vFun = LagrangianFunction.create(mesh, 2, obj.velocityOrder);
            obj.pFun = LagrangianFunction.create(mesh, 1, obj.pressureOrder);
        end

        function createPressureSmoother(obj)
            m              = obj.innerMesh;
            s.trial        = LagrangianFunction.create(m,1,obj.pressureOrder);
            s.filterType   = 'PDE';
            s.mesh         = m;
            s.boundaryType = 'Neumann';
            s.metric       = 'Isotropy';
            obj.filter     = Filter.create(s);
        end

        function createMaterial(obj)
            s.type       = 'STOKES';
            s.nelem      = obj.innerMesh.nelem;
            obj.material = Material.create(s);
        end

        function [isLeft,isBottom,isTop,isCyl] = computeBoundariesHandles(obj)
            xpos     = obj.center(1);
            ypos     = obj.center(2);
            r        = obj.radius;
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
            isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);
            isCyl    = @(coor) abs(abs(coor(:,1) - xpos).^2+abs(coor(:,2) - ypos).^2-r^2) < 5e-5;
        end

        function dirVel = computeDirichletConditionsHandles(obj,vInlet,vWalls)
            [isLeft,isBottom,isTop,isCyl] = obj.computeBoundariesHandles();

            dirVel{1}.domain    = @(coor) isLeft(coor) & not(isTop(coor) | isBottom(coor));
            dirVel{1}.direction = [1,2];
            dirVel{1}.value     = vInlet;

            dirVel{2}.domain    = @(coor) isTop(coor) | isBottom(coor);
            dirVel{2}.direction = [1,2];
            dirVel{2}.value     = vWalls;

            dirVel{3}.domain    = @(coor) isCyl(coor);
            dirVel{3}.direction = [1,2];
            dirVel{3}.value     = [0,0];
        end

        function [dirichlet,dir_dofs] = computeDirichletBCMatrices(obj,dirVel)
            velocityFun = obj.vFun;
            dirichlet   = [];
            dir_dofs    = [];
            for i = 1:length(dirVel)
                dirDofs = velocityFun.getDofsFromCondition(dirVel{i}.domain);
                nodes   = 1 + (dirDofs(2:2:end)-2)/velocityFun.ndimf;
                nodes2  = repmat(nodes, [1 2]);
                iNod    = sort(nodes2(:));
                mat12   = repmat([1;2], [length(iNod)/2 1]);
                valmat  = repmat(dirVel{i}.value', [length(iNod)/2 1]);

                dirichlet(size(dirichlet,1)+1:size(dirichlet,1)+length(iNod),:) = [iNod mat12 valmat];
                dir_dofs(size(dir_dofs,1)+1:size(dir_dofs,1)+length(iNod),1)    = dirDofs;
            end
        end

        function LHS = computeLHS(obj,dtime)
            mesh          = obj.innerMesh;
            s.type        = 'Stokes';
            s.dt          = dtime;
            s.mesh        = mesh;
            s.material    = obj.material;
            s.velocityFun = obj.vFun;
            s.pressureFun = obj.pFun;
            int           = LHSintegrator.create(s);
            LHS           = int.compute();
        end

        function forcesFormula = computeBodyForces(obj)
            s.fHandle     = @(coor) [0.*coor,0.*coor];
            s.ndimf       = 2;
            s.mesh        = obj.innerMesh;
            forcesFormula = AnalyticalFunction(s);
        end

        function RHS = computeRHS(obj,LHS,dirichlet,dir_dofs)
            forcesFormula   = obj.computeBodyForces();
            s.type          = 'Stokes';
            s.mesh          = obj.innerMesh;
            s.velocityFun   = obj.vFun;
            s.pressureFun   = obj.pFun;
            s.forcesFormula = forcesFormula;
            int             = RHSintegrator.create(s);
            F               = int.integrate();
            uD              = dirichlet(:,3);
            R               = -LHS(:,dir_dofs)*uD;
            RHS             = F + R;
        end

        function [vValues,pValues] = solveSystem(obj,LHS,RHS,dirichlet,dir_dofs)
            nDofs          = obj.vFun.nDofs + obj.pFun.nDofs;
            free_dofs_plus = setdiff(1:nDofs,dir_dofs);
            LHSr           = LHS(free_dofs_plus,free_dofs_plus);
            RHSr           = RHS(free_dofs_plus);
            x              = LHSr\RHSr;
            fullx          = obj.addDirichlet(x,dirichlet,dir_dofs,nDofs);
            ndofsV         = obj.vFun.nDofs;
            vValues        = fullx(1:ndofsV,:);
            pValues        = fullx(ndofsV+1:end,:);
        end

        function [v,p] = fillTrialFunctions(obj,vValues,pValues)
            v           = obj.vFun;
            pressureFun = obj.pFun;
            nu          = v.ndimf;
            nnode       = round(length(vValues)/nu);
            nodes       = 1:nnode;
            velfval     = zeros(nnode,nu);
            for idim = 1:nu
                dofs = nu*(nodes-1)+idim;
                velfval(:,idim) = vValues(dofs, end);
            end
            v.fValues = velfval;
            pressureFun.fValues = pValues(:,end);
            p = obj.filter.compute(pressureFun,'QUADRATIC');
        end
    end

    methods (Static, Access = private)
        function fullx = addDirichlet(x,dirichlet,dir_dofs,nDofs)
            uD                 = dirichlet(:,3);
            nsteps             = length(x(1,:));
            uD                 = repmat(uD,1,nsteps);
            fullx              = zeros(nDofs,nsteps);
            free_dofs          = setdiff(1:(nDofs),dir_dofs);
            fullx(free_dofs,:) = x;
            if ~isempty(dir_dofs)
                fullx(dir_dofs,:) = uD;
            end
        end
    end
end