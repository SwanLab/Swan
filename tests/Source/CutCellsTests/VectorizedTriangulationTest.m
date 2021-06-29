classdef VectorizedTriangulationTest < testShowingError
    
    properties (Access = protected)
        tol = 1e-14;
    end
    
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
        determinant
        determinant2
        
        connecBcutMesh
    end
    
    methods (Access = public)
        
        function obj = VectorizedTriangulationTest()
            obj.init();
            obj.computeDeterminant();
            obj.computeDeterminant2();
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
            obj.createUnfittedMesh();
            obj.plotUnfittedMeshCutPointsAndNormals();
            obj.computeValidInnerAndBoundaryCutMesh();            
        end
        
    end
    
    methods (Access = private)
        
        function createBackgroundMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            m = Mesh(s);
            obj.backgroundMesh = m;
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            errorV = obj.computeVolumeError();
            errorC = obj.computeConnecInnerCutMeshError();
            errorB = obj.computeConnecBoundaryCutMeshError();
            obj.error = errorV + errorC + errorB;
            det  = obj.determinant();
            det2 = obj.determinant2();
            bCutMeshVolume = obj.computeVolumes(obj.uMesh.boundaryCutMesh.mesh);
        end
        
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
            x = xCoord.xCutPoints(:,1);
            y = xCoord.xCutPoints(:,2);
            z = xCoord.xCutPoints(:,3);
            hold on
            plot3(x,y,z,'k*','LineWidth',10,'MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k')
        end
        
        function plotNormals(obj)
            m = obj.uMesh.boundaryCutMesh.mesh;            
            m.plotNormals();
        end
        
        function computeDeterminant2(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            mesh = Mesh(s);
            nDime   = mesh.ndim;
            nNode   = mesh.nnode;
            nElem   = mesh.nelem;
            q = Quadrature.set(mesh.type);
            q.computeQuadrature('LINEAR');
            mesh.interpolation.computeShapeDeriv(q.posgp)
            jacobian = zeros(nDime,nDime,nElem,q.ngaus);
            coordElem = permute(mesh.coordElem,[2 1 3]);
            for igaus = 1:q.ngaus
                dShapes = mesh.interpolation.deriv(:,:,igaus);
                jac = zeros(nDime,nDime,nElem);
                for kNode = 1:nNode
                    dShapeIK = dShapes(:,kNode);
                    xKJ      = coordElem(kNode,:,:);
                    jacIJ    = bsxfun(@times, dShapeIK, xKJ);
                    jac = jac + jacIJ;
                end
                jacobian(:,:,:,igaus) = jac;
            end
            executor = MatrixVectorizedInverterFactory().create(jacobian);
            obj.determinant2 = executor.computeDeterminant(jacobian);
        end
        
        function computeDeterminant(obj)
            xyz   = obj.coord;
            nodes = obj.connec;
            node1 = xyz(nodes(:,1),:);
            node2 = xyz(nodes(:,2),:);
            node3 = xyz(nodes(:,3),:);
            node4 = xyz(nodes(:,4),:);
            nelem = size(nodes,1);
            detA = zeros(nelem,1);
            for ielem = 1:nelem
                A = [node3(ielem,:) 1 ;node2(ielem,:) 1 ;node1(ielem,:) 1;node4(ielem,:) 1];
                detA(ielem) = det(A);
            end
            obj.determinant = detA;
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
            cU = obj.uMesh.boundaryCutMesh.mesh.connec;
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
            quad = Quadrature.set(mesh.type);
            quad.computeQuadrature('CONSTANT');
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
    
    methods (Access = protected, Abstract)
        init(obj)        
    end
    
end