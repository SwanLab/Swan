classdef StokesProblemBoundaryCondition < handle
    properties (Access = public)
        dirConditions
        dirDofs
        nodesConditions
    end

    properties (Access = private)
        height
        mesh
        uMesh
        velocityFun
        pressureFun
    end

    properties (Access = private)
        isLeft
        isRight
        isBottom
        isTop
        isHorizontal
        dirVelBC
        dirPreBC
        nnodesCutMesh
        dirDofsAirfoilMNodes
        middleNodesCoor
        dirDofsAirfoilNodes
        useMiddleNodes
    end

    methods (Access = public)
        function obj = StokesProblemBoundaryCondition(cParams)
            obj.init(cParams);
            obj.defineVariables();
            obj.defineDomainBoundaries();
            obj.useMiddleNodes = strcmpi(obj.velocityFun.order, 'P2');
        end

        function compute(obj)
            obj.setDomainDirichletVelocityBC();
            obj.setDofsAirfoil();
            obj.defineNodesConditions();
            obj.constructDirichletVelBCMatrix();
            obj.setDomainDirichletPressureBC();
            obj.constructDirichletPressureBCMatrix();
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.height      = cParams.height;
            obj.mesh        = cParams.mesh;
            obj.uMesh       = cParams.uMesh;
            obj.velocityFun = cParams.velocityFun;
            obj.pressureFun = cParams.pressureFun;
        end

        function defineVariables(obj)
            obj.dirConditions = [];
            obj.dirDofs       = [];
            obj.nodesConditions = [];
        end

        function defineDomainBoundaries(obj)
            tolerance = 1e-12;
            obj.isLeft       = @(coor) (abs(coor(:,1) - min(coor(:,1))) < tolerance);
            obj.isRight      = @(coor) (abs(coor(:,1) - max(coor(:,1))) < tolerance);
            obj.isBottom     = @(coor) (abs(coor(:,2) - min(coor(:,2))) < tolerance);
            obj.isTop        = @(coor) (abs(coor(:,2) - max(coor(:,2))) < tolerance);
            obj.isHorizontal = @(coor) (abs(coor(:,2) - obj.height/2) < tolerance);
        end

        function setDomainDirichletVelocityBC(obj)
            obj.dirVelBC{1}.domain    = @(coor) obj.isLeft(coor) & ~ (obj.isTop(coor) | obj.isBottom(coor));
            obj.dirVelBC{1}.direction = [1, 2];
            obj.dirVelBC{1}.value     = [1, 0];

            obj.dirVelBC{2}.domain    = @(coor) obj.isTop(coor) | obj.isBottom(coor);
            obj.dirVelBC{2}.direction = [1, 2];
            obj.dirVelBC{2}.value     = [0, 0];
        end

        function setDofsAirfoilNodes(obj)
            cutCoords = obj.uMesh.boundaryCutMesh.mesh.coord;
            tol = 1e-10;
            isAirfoil = @(coor) any(pdist2(coor, cutCoords) < tol, 2);
            obj.dirDofsAirfoilNodes = obj.velocityFun.getDofsFromCondition(isAirfoil);
            obj.dirDofsAirfoilNodes = sort(reshape(obj.dirDofsAirfoilNodes, [], 1));
        end

        function setAirfoilDirichletVelocityBC(obj)
            cutCoords = obj.uMesh.boundaryCutMesh.mesh.coord;
            ind = numel(obj.dirVelBC);
            obj.dirVelBC{ind+1}.domain = @(coor) any(ismember(coor, cutCoords, 'rows'), 2);
            obj.dirVelBC{ind+1}.direction = [1, 2];
            obj.dirVelBC{ind+1}.value = [0, 0];
        end

        function findMiddleNodes(obj)
            if ~obj.useMiddleNodes
                obj.middleNodesCoor = [];
                obj.dirVelBC{end+1}.value = [0, 0];
                return
            end
            bCutMesh = obj.uMesh.boundaryCutMesh.mesh;
            x1 = bCutMesh.coord(bCutMesh.connec(:,1), :);
            x2 = bCutMesh.coord(bCutMesh.connec(:,2), :);
            obj.middleNodesCoor = (x1 + x2) / 2;
        end

        function setDofsAirfoilMNodes(obj)
            if ~obj.useMiddleNodes
                obj.dirDofsAirfoilMNodes = [];
                return
            end
            nEdges = size(obj.uMesh.boundaryCutMesh.mesh.connec, 1);
            obj.dirDofsAirfoilMNodes = zeros(2, nEdges);
            for i = 1:nEdges
                isxcoord = @(coor) abs(coor(:,1) - obj.middleNodesCoor(i,1)) < 1e-12;
                isycoord = @(coor) abs(coor(:,2) - obj.middleNodesCoor(i,2)) < 1e-12;
                isAirfoil = @(coor) isxcoord(coor) & isycoord(coor);

                obj.dirDofsAirfoilMNodes(:, i) = obj.velocityFun.getDofsFromCondition(isAirfoil);
            end
            obj.dirDofsAirfoilMNodes = sort(reshape(obj.dirDofsAirfoilMNodes, [], 1));
            obj.dirVelBC{end+1}.value = [0, 0];
        end

        function setDofsAirfoil(obj)
            obj.setDofsAirfoilNodes();
            obj.setAirfoilDirichletVelocityBC();
            obj.findMiddleNodes();
            obj.setDofsAirfoilMNodes();
        end

        function defineNodesConditions(obj)
            obj.nodesConditions = 1 + (obj.dirDofsAirfoilNodes(2:2:end)-2)/obj.velocityFun.ndimf;
        end

        function constructDirichletVelBCMatrix(obj)
            n = numel(obj.dirVelBC);
            for i = 1:n
                if i < n-1
                    dirDofsVel = obj.velocityFun.getDofsFromCondition(obj.dirVelBC{i}.domain);
                elseif i == n-1
                    dirDofsVel = obj.dirDofsAirfoilNodes;
                else
                    if ~isempty(obj.dirDofsAirfoilMNodes)
                        dirDofsVel = obj.dirDofsAirfoilMNodes;
                    else

                        continue;
                    end
                end

                nodes = 1 + (dirDofsVel(2:2:end)-2)/obj.velocityFun.ndimf;
                nodesDuplicated = repmat(nodes, [1 2]);
                sortedNodes = sort(nodesDuplicated(:));
                direction = repmat([1; 2], [length(sortedNodes)/2 1]);
                value = repmat(obj.dirVelBC{i}.value', [length(sortedNodes)/2 1]);

                obj.dirConditions = [obj.dirConditions; [sortedNodes direction value]];
                obj.dirDofs = [obj.dirDofs; dirDofsVel];
            end
        end


        function setDomainDirichletPressureBC(obj)
            obj.dirPreBC{1}.domain = @(coor) obj.isRight(coor);
            obj.dirPreBC{1}.direction = 1;
            obj.dirPreBC{1}.value = 0;
        end

        function constructDirichletPressureBCMatrix(obj)
            for i = 1:length(obj.dirPreBC)
                dirDofsPre = obj.pressureFun.getDofsFromCondition(obj.dirPreBC{i}.domain);
                dirDofsPre = dirDofsPre(ceil(length(dirDofsPre)/2));
                direction = ones(size(dirDofsPre));
                value = ones(size(dirDofsPre)) * obj.dirPreBC{i}.value;
                obj.dirConditions(end+1:end+length(dirDofsPre), :) = [dirDofsPre + obj.velocityFun.nDofs direction value];
                obj.dirDofs(end+1:end+length(dirDofsPre), 1) = dirDofsPre + obj.velocityFun.nDofs;
            end
        end

    end
end
