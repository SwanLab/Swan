classdef TestNaca < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        refMesh
        levelSet
        rawMesh
        mesh
    end
    
    methods (Access = public)
        
        function obj = TestNaca()
            close all;
            obj.createReferenceMesh();
            obj.levelSet = obj.createLevelSet(obj.refMesh);
            obj.createFluidMesh();
            obj.createFluidMeshGoodConditioning();
        end
        
    end
    
    methods (Access = private)
                
        function createReferenceMesh(obj)      
            obj.refMesh = TriangleMesh(2,1,150,75);
        end
        
        function createFluidMesh(obj)
            s.backgroundMesh = obj.refMesh;
            s.boundaryMesh   = obj.refMesh.createBoundaryMesh();
            uMesh            = UnfittedMesh(s);
            uMesh.compute(obj.levelSet.fValues);       
            obj.rawMesh = uMesh.createInnerMesh();
        end

        function createFluidMeshGoodConditioning(obj)
            m        = obj.createMeshAlphaTriangulation();
            s.connec = obj.computeConnectivitiesGoodCond(m);
            s.coord  = obj.rawMesh.coord;
            m2       = Mesh.create(s);
            obj.mesh = m2.computeCanonicalMesh();
            figure;
            obj.mesh.plot();
        end

        function m = createMeshAlphaTriangulation(obj)
            points = obj.rawMesh.coord;
            T      = alphaShape(points,.007);
            DT     = alphaTriangulation(T);

            s.connec = DT;
            s.coord  = points;
            m        = Mesh.create(s);
        end

        function connec = computeConnectivitiesGoodCond(obj,mAlpha)
            xV(:,1,:) = mAlpha.computeBaricenter();
            lsFun  = obj.createLevelSet(mAlpha);
            lsElem = squeeze(lsFun.evaluate(xV));
            connec = mAlpha.connec(lsElem<=0,:);
        end

    end

    methods (Static, Access = private)

        function ls = createLevelSet(m)
            s.type = 'Naca';
            s.xLE  = 0.5;
            s.yLE  = 0.5;

            s.chord = 1;
            s.p     = 0.5;
            s.m     = 0.02;
            s.t     = 0.12;

            g  = GeometricalFunction(s);
            ls = g.computeLevelSetFunction(m);
        end

    end

end