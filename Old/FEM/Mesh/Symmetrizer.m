classdef Symmetrizer < handle
    
    properties (Access = public)
        symmetricMesh
    end
    
    properties (Access = private)
        symmetricPoint
        meshMerger
    end
    
    properties (Access = private)
        mesh
        symmetrizedMesh
        symmetricLine
    end
    
    methods (Access = public)
        
        function obj = Symmetrizer(cParams)
            obj.init(cParams);
            obj.createSymmetricPointComputer();
            obj.createSymmetrizedMesh();
            obj.createMeshMerger();
        end
        
        function m = computeSymmetricMesh(obj)
            m = obj.meshMerger.compute();
        end
        
        function sF = symmetrizeScalarField(obj,f)
            nodes = obj.meshMerger.computeRemainingNodes();
            fS = f;
            fD = [f;fS];
            sF   = fD(nodes,:);
        end
        
        function sF = symmetrizeVectorField(obj,f)
            nodes = obj.meshMerger.computeRemainingNodes();
            fS = obj.computeFsymmetric(f);
            fD   = [f;fS];
            sF   = fD(nodes,:);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.symmetricLine = cParams.symmetricLine;
        end
        
        function fSym = computeFsymmetric(obj,f)
            s.vector = f;
            s.line  = obj.symmetricLine;
            s.line.point = [0,0];
            sP = SymmetricPointComputer(s);
            fSym = sP.computeSymmetricVector();
        end
        
        function createSymmetricPointComputer(obj)
            s.vector = obj.mesh.coord;
            s.line   = obj.symmetricLine;
            sP = SymmetricPointComputer(s);
            obj.symmetricPoint = sP;
        end
        
        function createSymmetrizedMesh(obj)
            sConnec(:,[1 3 2]) = obj.mesh.connec(:,[1 2 3]);
            sCoord = obj.symmetricPoint.computeSymmetricVector();
            s.connec = sConnec;
            s.coord  = sCoord;
            m = Mesh.create(s);
            obj.symmetrizedMesh = m;
        end

        function createMeshMerger(obj)
            isInSymLine = obj.symmetricPoint.isNodeInLine;
            s.meshA  = obj.mesh;
            s.meshB  = obj.symmetrizedMesh;
            s.isMergedNodeA = isInSymLine;
            s.isMergedNodeB = isInSymLine;
            obj.meshMerger = MeshMerger(s);
        end
        
    end
    
end