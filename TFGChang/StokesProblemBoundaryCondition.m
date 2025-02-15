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
        lengthCutMesh
        dirDofsAirfoilMNodes
        middleNodesCoor
        dirDofsAirfoilNodes
    end

    methods (Access = public)
        
        function obj = StokesProblemBoundaryCondition(cParams)
            obj.init(cParams);
            obj.defineVariables();
            obj.defineDomainBoundaries();
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

        function init(obj,cParams)
            obj.height         = cParams.height;
            obj.mesh           = cParams.mesh;
            obj.uMesh          = cParams.uMesh;
            obj.velocityFun    = cParams.velocityFun;
            obj.pressureFun    = cParams.pressureFun;
        end

        function defineVariables(obj)
            obj.lengthCutMesh         = size(obj.uMesh.boundaryCutMesh.mesh.coord,1);
            obj.dirDofsAirfoilNodes   = zeros(2,obj.lengthCutMesh);
            obj.middleNodesCoor       = zeros(obj.lengthCutMesh,2);
            obj.dirDofsAirfoilMNodes  = zeros(2,obj.lengthCutMesh);
            obj.dirConditions         = [];
            obj.dirDofs               = [];
        end  

        function defineDomainBoundaries(obj)
            tolerance        = 1e-12;
            obj.isLeft       = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < tolerance);
            obj.isRight      = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < tolerance);
            obj.isBottom     = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < tolerance);
            obj.isTop        = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < tolerance);
            obj.isHorizontal = @(coor) (abs(coor(:,2) - obj.height/2) < tolerance); 
        end

        function setDomainDirichletVelocityBC(obj)
            obj.dirVelBC{2}.domain    = @(coor) obj.isTop(coor) | obj.isBottom(coor);
            obj.dirVelBC{2}.direction = [1,2];
            obj.dirVelBC{2}.value     = [0,0]; 
            obj.dirVelBC{1}.domain    = @(coor) obj.isLeft(coor) & not(obj.isTop(coor) | obj.isBottom(coor));
            obj.dirVelBC{1}.direction = [1,2];
            obj.dirVelBC{1}.value     = [1,0]; 
        end

        function setDofsAirfoilNodes(obj)
            cutCoords = obj.uMesh.boundaryCutMesh.mesh.coord(1:obj.lengthCutMesh, :);
            isAirfoil = @(coor) any(ismember(coor, cutCoords, 'rows'), 2);
            obj.dirDofsAirfoilNodes = obj.velocityFun.getDofsFromCondition(isAirfoil);
            obj.dirDofsAirfoilNodes = sort(reshape(obj.dirDofsAirfoilNodes, [], 1));
        end

        function findMiddleNodes(obj)
            obj.middleNodesCoor =(obj.uMesh.boundaryCutMesh.mesh.coord(obj.uMesh.boundaryCutMesh.mesh.connec(:,1),:) + obj.uMesh.boundaryCutMesh.mesh.coord(obj.uMesh.boundaryCutMesh.mesh.connec(:,2),:))./2;
        end

        function setDofsAirfoilMNodes(obj)
            for i = 1:1:obj.lengthCutMesh
                isxcoord    = @(coor) coor(:,1) == obj.middleNodesCoor(i,1);
                isycoord    = @(coor) coor(:,2) == obj.middleNodesCoor(i,2);
                isAirfoil  = @(coor) isxcoord(coor) & isycoord(coor);
            
                obj.dirDofsAirfoilMNodes(:,i) = obj.velocityFun.getDofsFromCondition(isAirfoil);
            end
            obj.dirDofsAirfoilMNodes = sort(reshape(obj.dirDofsAirfoilMNodes,size(obj.dirDofsAirfoilMNodes,2)*2,1));
        end

        function setDofsAirfoil(obj)
            obj.setDofsAirfoilNodes();
            obj.findMiddleNodes();
            obj.setDofsAirfoilMNodes();
        end

        function defineNodesConditions(obj)
            obj.nodesConditions = 1 + (obj.dirDofsAirfoilNodes(2:2:end)-2)/obj.velocityFun.ndimf;
        end

        function constructDirichletVelBCMatrix(obj)
            for i = 1:1:4 
                if i == 1 || i == 2
                    dirDofsVel = obj.velocityFun.getDofsFromCondition(obj.dirVelBC{i}.domain);
                elseif i == 3
                    dirDofsVel = obj.dirDofsAirfoilNodes;
                    obj.dirVelBC{i}.value     = [0,0]; 
                elseif i == 4
                    dirDofsVel = obj.dirDofsAirfoilMNodes;
                    obj.dirVelBC{i}.value     = [0,0]; 
                end
                
                nodes            = 1 + (dirDofsVel(2:2:end)-2)/obj.velocityFun.ndimf;
                nodesDuplicated  = repmat(nodes, [1 2]);
                sortedNodes      = sort(nodesDuplicated(:));
                direction        = repmat([1;2], [length(sortedNodes)/2 1]);
                value            = repmat(obj.dirVelBC{i}.value', [length(sortedNodes)/2 1]);
                obj.dirConditions(size(obj.dirConditions,1)+1:size(obj.dirConditions,1)+length(sortedNodes),:) = [sortedNodes direction value];
                obj.dirDofs(size(obj.dirDofs,1)+1:size(obj.dirDofs,1)+length(sortedNodes),1)                   = dirDofsVel;
            end
        end

        function setDomainDirichletPressureBC(obj)
            obj.dirPreBC{1}.domain    = @(coor) obj.isRight(coor) & obj.isHorizontal(coor);
            obj.dirPreBC{1}.direction = 1;
            obj.dirPreBC{1}.value     = 0;
        end

        function constructDirichletPressureBCMatrix(obj)
            for i = 1:length(obj.dirPreBC)
                dirDofsPre  = obj.pressureFun.getDofsFromCondition(obj.dirPreBC{i}.domain);
                direction   = ones(size(dirDofsPre));
                value       = ones(size(dirDofsPre)).*obj.dirPreBC{i}.value';
                obj.dirConditions(size(obj.dirConditions,1)+1:size(obj.dirConditions,1)+length(dirDofsPre),:) = [dirDofsPre+obj.velocityFun.nDofs direction value];
                obj.dirDofs(size(obj.dirDofs,1)+1:size(obj.dirDofs,1)+length(dirDofsPre),1)                   = dirDofsPre+obj.velocityFun.nDofs;
            end
        end
        
    end

end