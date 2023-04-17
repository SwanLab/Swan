classdef LevelSetYounes < LevelSetCreator

    properties (Access = private)
        epsilon
        cellCoord
        phi
        y1
        y2
        deformedCoord
        m1
        m2
        mesh
        remesher
        cellLevelSetParams
        orientationVectors
        nCells
    end
    
    methods (Access = public)
        
        function obj = LevelSetYounes(cParams)
            obj.init(cParams);
            obj.createDeformedCoord();
            obj.createRemesher();
            obj.interpolateDeformedCoord();
            obj.createM1M2();
            obj.interpolateM1M2();
        end

        function ls = computeLS(obj,epsilons)
            nEps = length(epsilons);
            ls = cell(nEps,1);
            for iEps = 1:nEps
                obj.epsilon = epsilons(iEps);
                obj.computeLevelSet();
                ls{iEps} = obj.getValue();
            end
        end

        function mF = getFineMesh(obj)
            fMesh = obj.remesher.fineMesh;
            mF = fMesh.createDiscontinuousMesh();
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.createCellCoord();
            % PAO: thresholdParameters
            obj.createCellLevelSet();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.orientationVectors = cParams.orientationVectors;
            obj.cellLevelSetParams = cParams.cellLevelSetParams;
        end

        function createDeformedCoord(obj)
            c  = obj.orientationVectors;
            dC = c.computeDeformedCoordinates();
            obj.deformedCoord = dC;
        end

        function createRemesher(obj)
            s.mesh    = obj.mesh.createDiscontinuousMesh();
            s.nLevels = 2;
            r  = Remesher(s);
            r.remesh();
            obj.remesher = r;
         end

        function interpolateDeformedCoord(obj)
            y = obj.deformedCoord;
            y = obj.interpolateDiscontinousFunction(y);
            y = abs(y);
            obj.deformedCoord.fValues = y;
        end

        function createM1M2(obj)
            x1 = obj.mesh.coord(:,1);
            x2 = obj.mesh.coord(:,2);
            obj.m1 = 1-x1.*x2;
            obj.m2 = 1-x1.*x2;
        end
        
        function fDI = interpolateContinousFunctionToDisc(obj,fC)
            fD   = obj.createDiscontinousValues(fC);
            fDI  = obj.interpolateDiscontinousFunction(fD);
        end    

        function interpolateM1M2(obj)
            m1 = obj.m1;
            m2 = obj.m2;
            obj.m1 = abs(obj.interpolateContinousFunctionToDisc(m1));
            obj.m2 = abs(obj.interpolateContinousFunctionToDisc(m2));
        end

        function fV = createDiscontinousValues(obj,f)
            s.mesh    = obj.mesh;
            s.fValues = f;
            fC = P1Function(s);
            fD = fC.project('P1D');
            fV = fD;
        end

        function vq = interpolateDiscontinousFunction(obj,v)
            f = v;
            r = obj.remesher;
            f = r.interpolate(f);
            vq = f.getFvaluesAsVector();
        end 

        function createCellCoord(obj)
            nDim = obj.mesh.ndim;
            for iDim = 1:nDim
                x = obj.deformedCoord.fValues(:,iDim);
                y = obj.computeMicroCoordinate(x);
                y = obj.periodicFunction(y);
                yT(:,iDim) = y;
            end
            obj.cellCoord = yT;
        end

         function y = computeMicroCoordinate(obj,x)
            eps = obj.epsilon;
            y = (x-min(x))/eps;
         end
        
        function createCellLevelSet(obj)
            %s       = obj.cellLevelSetParams;
            s.coord = obj.cellCoord;
            obj.computeLevelSetValue(s);
        end

        function computeLevelSetValue(obj,cParams)
            q=2;

            fpy1 = cParams.coord(:,1);
            fpy2 = cParams.coord(:,1);

            ls(:,1) = ((2*fpy1)./obj.m1).^q+((2*fpy2)./obj.m2).^q;
            %ls(:,2) = ls(:,1);
            obj.levelSet = 1-ls;
        end
    end

    methods (Access = private, Static)

        function f = periodicFunction(y)
            %f = abs(cos(pi/2*(y))).^2;
            f = y - floor(y) - 1/2;
        end

    end
        
    
end
