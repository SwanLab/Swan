classdef AeroForcesCalculation < handle
    
    properties (Access = public)
       L
       D
       E
    end
    
    properties (Access = private)
        mesh
        nodesConditions
        pressureFun
    end

    properties (Access = private)
        connec
        bMesh
        pBoundary
        normalVectors
        lengthElement
        centroid
        centralPointsElement
        nX
        nY
    end

    methods (Access = public)
        
        function obj = AeroForcesCalculation(cParams)
            obj.init(cParams);
            obj.identifyBoundaryEdges();
            obj.createBoundaryMesh();
            obj.assignPressureToBoundary();
            obj.defineVariables();
        end

        function compute(obj)
            obj.plotPressureDistribution();
            obj.computeNormalVectors();
            obj.convertNorVectorsToLanFunctions();
            obj.computeDragForces();
            obj.computeLiftForces();
            obj.computeEfficiency();
            obj.plotResults();
        end
        
    end
        
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.nodesConditions   = cParams.nodesConditions;
            obj.pressureFun       = cParams.pressureFun;
        end

        function identifyBoundaryEdges(obj)
            obj.mesh.computeEdges();
            e  = obj.mesh.edges.nodesInEdges;
            bE = ismember(e,obj.nodesConditions);
            bEIndices = find(prod(bE,2));
            obj.connec = e(bEIndices,:);
        end

        function createBoundaryMesh(obj)
            ss.coord    = obj.mesh.coord;
            ss.connec   = obj.connec;
            ss.kFace    = -1;
            obj.bMesh       = Mesh.create(ss);
            obj.bMesh       = obj.bMesh.computeCanonicalMesh();
        end

        function assignPressureToBoundary(obj)
            obj.pBoundary         = LagrangianFunction.create(obj.bMesh,1,obj.pressureFun.order); 
            obj.pBoundary.fValues = obj.pressureFun.fValues(obj.nodesConditions,1);
        end

        function plotPressureDistribution(obj)
            obj.pBoundary.plot();
        end

        function defineVariables(obj)
            obj.normalVectors        = zeros(obj.bMesh.nelem,obj.bMesh.ndim);
            obj.lengthElement        = zeros(obj.bMesh.nelem,1);
            obj.centroid             = mean(obj.bMesh.coord);
            obj.centralPointsElement = (obj.bMesh.coord(obj.bMesh.connec(:,1),:) + obj.bMesh.coord(obj.bMesh.connec(:,2),:))/2;
        end

        function computeNormalVectors(obj)
            cont = 1;

            for iE = 1:obj.bMesh.nelem
                node1 = obj.bMesh.coord(obj.bMesh.connec(iE,1),:);
                node2 = obj.bMesh.coord(obj.bMesh.connec(iE,2),:);
            
                if node1(1)<= 5
                tanVect                    = (node2-node1)/(abs(norm(node2-node1))); 
                obj.normalVectors(cont,:) = -tanVect * [0 -1;1 0];
                obj.lengthElement(cont)    = abs(norm(node1-node2));
            
                cont = cont +1;
                end
            
            end
        end

        function convertNorVectorsToLanFunctions(obj)
            obj.nX = LagrangianFunction.create(obj.bMesh,1,'P0');
            obj.nY = LagrangianFunction.create(obj.bMesh,1,'P0');
            obj.nX.fValues = obj.normalVectors(:,1);
            obj.nY.fValues = obj.normalVectors(:,2);
        end
        
        function computeDragForces(obj)
            sss.operation    = @(x) -obj.pBoundary.evaluate(x).*obj.nX.evaluate(x);
            pNX              = DomainFunction(sss);
            obj.D            = Integrator.compute(pNX,obj.bMesh,2);
        end

        function computeLiftForces(obj)
            sss.operation = @(x) -obj.pBoundary.evaluate(x).*obj.nY.evaluate(x);
            pNY           = DomainFunction(sss);
            obj.L         = Integrator.compute(pNY,obj.bMesh,2);
        end

        function computeEfficiency(obj)
            obj.E = obj.L/obj.D;
        end

        function plotResults(obj)
            quiver(obj.centralPointsElement(:,1),obj.centralPointsElement(:,2),obj.normalVectors(:,1),obj.normalVectors(:,2));
            hold on;
            quiver(obj.centroid(1,1),obj.centroid(1,2),obj.D,0);
            hold on;
            quiver(obj.centroid(1,1),obj.centroid(1,2),0,obj.L);
            hold on;
            obj.bMesh.plot() 
        end
    end

end