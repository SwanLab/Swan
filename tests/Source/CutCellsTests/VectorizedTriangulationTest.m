classdef VectorizedTriangulationTest < handle
    
    properties (Access = private)
        uMesh
        validInnerCutMesh
        backgroundMesh
        boundaryMesh
    end
    
    properties (Access = protected)
        coord
        connec
        boundaryConnec
        levelSet
        connecBcutMesh
    end
    
    methods (Access = public)
        
        function obj = VectorizedTriangulationTest(cParams)
            obj.init(cParams);
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
            obj.createUnfittedMesh();
            obj.plotUnfittedMeshCutPointsAndNormals();
            obj.computeValidInnerAndBoundaryCutMesh();
        end
        
        function error = computeError(obj)
            errorV = obj.computeVolumeError();
            errorC = obj.computeConnecInnerCutMeshError();
            errorB = obj.computeConnecBoundaryCutMeshError();
            error = errorV + errorC + errorB;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.coord    = cParams.coord;
            obj.connec   = cParams.connec;
            obj.levelSet = cParams.levelSet;
            obj.boundaryConnec = cParams.boundaryConnec;
            obj.connecBcutMesh = cParams.connecBcutMesh;
        end
        
        function createBackgroundMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            m = Mesh.create(s);
            obj.backgroundMesh = m;
        end

    end
    
    methods (Access = protected)
        
        function computeValidInnerAndBoundaryCutMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.levelSet       = obj.levelSet;
            s.boundaryConnec = obj.boundaryConnec;
            c = ComputingInnerAndBoundaryCutMesh(s);
            c.compute();
            obj.validInnerCutMesh = c.innerCutMesh;
        end
        
        function createBoundaryMesh(obj)
            bConnec = obj.boundaryConnec;
            for iFace = 1:size(bConnec,1)
                con = bConnec(iFace,:);
                s.coord =  obj.backgroundMesh.coord(con,:);
                s.nodesInBoxFaces = con;
                s.connec = [1 2 3];
                s.kFace = 0;
                s.dimension = [];
                s.isRectangularBox = false;
                m{iFace} = BoundaryMesh(s);
            end
            obj.boundaryMesh = m;
        end
        
        function createUnfittedMesh(obj)
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.backgroundMesh;
            obj.uMesh = UnfittedMesh(s);
            obj.uMesh.compute(obj.levelSet)
        end
        
        function plotUnfittedMeshCutPointsAndNormals(obj)
            figure()
            obj.uMesh.plotBoundary();
            view([1 1 1])
            obj.plotCutPoints();
            camlight
            lighting flat
            obj.plotNormals();
        end
        
        
        function plotCutPoints(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.levelSet       = obj.levelSet;
            s.boundaryConnec = obj.boundaryConnec;
            c = ComputingInnerAndBoundaryCutMesh(s);
            c.compute();
            xCoord = c.cutCoordComputer;
            x = xCoord.cutValues(:,1);
            y = xCoord.cutValues(:,2);
            z = xCoord.cutValues(:,3);
            hold on
            plot3(x,y,z,'k*','LineWidth',10,'MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k')
        end
        
        function plotNormals(obj)
            m = obj.uMesh.boundaryCutMesh.mesh;
            m.plotNormals();
        end
        
        function error = computeConnecInnerCutMeshError(obj)
            cV = obj.validInnerCutMesh.mesh.connec;
            cU = obj.uMesh.innerCutMesh.mesh.connec;
            
            cV = obj.sortConnec(cV);
            cU = obj.sortConnec(cU);
            error = norm(cV(:)-cU(:));
        end
        
        function error = computeConnecBoundaryCutMeshError(obj)
            cV = obj.connecBcutMesh;
            cU = obj.uMesh.boundaryCutMesh.mesh.connec; % peta aqui
            areEquiv = obj.areConnecEquivalent(cU,cV);
            error = sum(~areEquiv)/length(areEquiv);
        end
        
        function error = computeVolumeError(obj)
            vV = obj.computeVolumes(obj.validInnerCutMesh.mesh);
            vU = obj.computeVolumes(obj.uMesh.innerCutMesh.mesh);
            vAT = zeros(size(vV));
            vAT(1:length(vU)) = vU;
            volums = [vV; vAT]';
            error  = abs(sum(vU) - sum(vV));
        end
        
    end
    
    methods (Access = private, Static)
        
        function v = computeVolumes(mesh)
            quad = Quadrature.create(mesh,0);
            v = mesh.computeDvolume(quad);
        end
        
        function c = sortConnec(connec)
            c = sort(connec')';
        end
        
        function areEquiv = areConnecEquivalent(connecA,connecB)
            isEqual = false(size(connecA));
            for iperm = 1:size(connecA,2)
               connecBp = circshift(connecB',iperm)';
               isEqual(:,iperm) = sum(abs(connecA - connecBp),2) < 1e-14;
            end
            areEquiv = any(isEqual,2);
        end
        
    end

end

