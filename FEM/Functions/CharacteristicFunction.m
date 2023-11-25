classdef CharacteristicFunction < L2Function

    properties (Access = private)
        coorP1
    end

    properties (Access = public)
        ndimf
    end

    properties (Access = private)
        fxy
        levelSet
        unfittedMesh
    end

    methods (Access = public)

        function obj = CharacteristicFunction(cParams)
            obj.init(cParams);
            obj.createP1CoorFunction();
            obj.createLevelSetFunction();
            obj.createUnfittedMesh();
%             obj.createRHS();
        end

        function fxV = evaluate(obj,xV)
            xy    = obj.coorP1.evaluate(xV);
            nGaus = size(xy,2);
            nElem = size(xy,3);
            fxV   = zeros(1,nGaus,nElem);
            for iGaus = 1:nGaus
                x = xy(1,iGaus,:);
                y = xy(2,iGaus,:);
                if (size(xy,1) == 3)
                    z = xy(3,iGaus,:);
                    f = obj.fxy(x,y,z);
                else
                    f = obj.fxy(x,y);
                end
                fxV(1,iGaus,f>0) = 0;
                fxV(1,iGaus,f<=0) = 1;
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ndimf = 1;
            obj.mesh  = cParams.mesh;
            obj.fxy   = cParams.fxy;
        end

        function createP1CoorFunction(obj)
            s.mesh     = obj.mesh;
            s.fValues  = obj.mesh.coord;
            obj.coorP1 = P1Function(s);
        end

        function createLevelSetFunction(obj)
%             fxy = @(x,y) (x-0.5).^2+(y-0.5).^2-0.3.^2;
            fxy = @(x,y,z) -((x-0.5).^2+(y-0.5).^2+(z-0.5).^2 -0.3.^2);
            xy    = obj.coorP1.fValues;
            val = fxy(xy(:,1), xy(:,2), xy(:,3));
            a.fValues = val;
            a.mesh = obj.mesh;
            lS = P1Function(a);
            lS.plot
            obj.levelSet = lS;
        end

        function createUnfittedMesh(obj)
%             b.levelSet = obj.levelSet.fValues;
            b.backgroundMesh = obj.mesh;
            b.boundaryMesh = obj.mesh.createBoundaryMesh();
            uM = UnfittedMesh(b);
            uM.compute(obj.levelSet.fValues);
            obj.unfittedMesh = uM;

            h = obj.unfittedMesh.innerMesh.mesh.computeMeanCellSize();
            s.type = 'Matlab'; % GiD
            s.filename = 'provant123';
            s.meshElementSize = num2str(h);
            s.meshFileName    = 'hmmmm22';
            s.swanPath        = '/home/ton/Github/Swan/';
            s.gidPath         = '/home/ton/GiDx64/gid-16.1.2d/';

            innermesh = uM.createFullInnerMesh(s);
        end
% 
%         function createRHS(obj)
%             uMesh = obj.unfittedMesh;
%             s.mesh = uMesh;
%             s.type = 'Unfitted';
%             int = RHSintegrator.create(s);
%             ls = obj.levelSet.fValues;
%             fNodes = ones(size(ls));
%             fInt = int.integrateInDomain(fNodes);
%         end

    end
end