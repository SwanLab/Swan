classdef Tryinggg < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        meshDomain
        nSubdomains
        tCross
        tFrame
        xmax
        xmin
        ymax
        ymin
        cellMesh
        referenceMesh
        nelem
        tolSameNode
    end

    methods (Access = public)

        function obj = Tryinggg()
            obj.init()

        end

    end

    methods (Access = private)

        function init(obj)
            s.Training  = [];                 % 'EIFEM'/'Multiscale'/'EIFisol'
            s.Option    = 'Dataset';          % 'Dataset'/'NN'
            obj.nelem     =  30;                %  Mesh refining

            obj.tolSameNode = 1e-11;


            %NON-UNIFORM DISTRIBUTION
            obj.tFrame= [0.15, 0.30, 0.45, 0.10, 0.25, 0.40, 0.35, 0.20, 0.15, 0.45;
                0.40, 0.10, 0.25, 0.35, 0.15, 0.45, 0.20, 0.30, 0.40, 0.10;
                0.25, 0.45, 0.15, 0.20, 0.40, 0.30, 0.10, 0.35, 0.25, 0.35];

            obj.tCross=[0.10, 0.45, 0.60, 0.25, 0.35, 0.50, 0.15, 0.40, 0.20, 0.55;
                0.50, 0.20, 0.30, 0.60, 0.10, 0.45, 0.55, 0.15, 0.35, 0.25;
                0.35, 0.55, 0.15, 0.40, 0.60, 0.20, 0.45, 0.30, 0.50, 0.10];

            % MULTISCALE ISOLATED
            s.Training  = 'Multiscale';
            s.Sampling  = 'Isolated';

            obj.nSubdomains = size(obj.tFrame');

            mSbd = obj.createSubDomainMeshes();
            [mD,mSb,iC,lG,iCR,discmesh] = obj.createMeshDomainJoiner(mSbd);
            obj.meshDomain = mD;        % mD:conj subdominis --> Tot el domini

            ls= obj.computeLevelSet();
            sUm.backgroundMesh = obj.meshDomain;
            sUm.boundaryMesh   = obj.meshDomain.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(-ls);

      %      uMesh.plot()
            aaa=uMesh.createInnerMesh();
            aaa.plot()


        end


        function  mSbd = createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = 2; Ly = 2;
            cM=cell(nY,nX);
            mSbd=cell(nY,nX);
            for jDom = 1:nY
                for iDom = 1:nX
                    refMesh  =obj.createStructuredMesh();
                    coord0  = refMesh.coord;
                    s.coord(:,1) = coord0(:,1)+Lx*(iDom-1);
                    s.coord(:,2) = coord0(:,2)+Ly*(jDom-1);
                    s.connec = refMesh.connec;
                    mIJ     = Mesh.create(s);
                    %                     plot(mIJ)
                    %                     hold on;
                    mSbd{jDom,iDom} = mIJ;
                    cM{jDom,iDom} = refMesh;  %same but with local coordinates
                end
            end
            obj.referenceMesh = mSbd{1,1};
            obj.cellMesh=cM;
        end

        function mS = createStructuredMesh(obj)
            n =obj.nelem;
            x1      = linspace(-1,1,n);
            x2      = linspace(-1,1,n);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            %
            % mesh= QuadMesh(1,1,n,n);
            % s.coord= mesh.coord;
            % s.connec=mesh.connec;

            obj.xmin = min(x1);
            obj.xmax = max(x1);
            obj.ymin = min(x2);
            obj.ymax = max(x2);
            delta = 1e-8;
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[delta,-delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[delta,delta];
            mS = Mesh.create(s);
        end


        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomainJoiner(obj,mSbd)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.referenceMesh;
            s.tolSameNode   = obj.tolSameNode;
            s.meshSbd       = mSbd;
            m = MeshJoiner(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        function ls=computeLevelSet(obj)
            [x0,y0] = obj.computeSubdomainCentroid();
            [Nx,Ny] = size(obj.tFrame);
            GeomParams(Nx,Ny) = struct('type',[],'length',[],'tFrame',[],'tCross',[],'xCoorCenter',[],'yCoorCenter',[]);

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    GeomParams(i,j).type        = "CrossedSquare";
                    GeomParams(i,j).length      = 2;
                    GeomParams(i,j).tFrame      = obj.tFrame(i,j);
                    GeomParams(i,j).tCross      = obj.tCross(i,j);
                    GeomParams(i,j).xCoorCenter = x0(i,j);
                    GeomParams(i,j).yCoorCenter = y0(i,j);

                    % geomFun = GeometricalFunction(GeomParams(1,1));
                    % ls = geomFun.computeLevelSetFunction(obj.meshDomain);
                    % 
                    % sUm.backgroundMesh = obj.meshDomain;
                    % sUm.boundaryMesh   = obj.meshDomain.createBoundaryMesh;
                    % uMesh              = UnfittedMesh(sUm);
                    % uMesh.compute(-ls.fValues);
                    % 
                    % %      uMesh.plot()
                    % aaa=uMesh.createInnerMesh();
                    % aaa.plot()
                end
            end
            s.type        = 'GivenPattern';
            s.paramsList  = GeomParams;
            g             = GeometricalFunction(s);
            phiFun        = g.computeLevelSetFunction(obj.meshDomain);
            ls            = phiFun.fValues;
        end

        function [x0,y0]= computeSubdomainCentroid(obj)
            % [Nx, Ny] = size(obj.r);
            xMin=min(obj.meshDomain.coord(:,1));
            xMax=max(obj.meshDomain.coord(:,1));
            yMin=min(obj.meshDomain.coord(:,2));
            yMax=max(obj.meshDomain.coord(:,2));

            dx = obj.xmax-obj.xmin;
            dy = obj.ymax-obj.ymin;

            x_center = xMin + dx/2 : dx : xMax - dx/2;
            y_center = yMin + dy/2 : dy : yMax - dy/2;

            [x0, y0] = meshgrid(x_center, y_center);
        end
    end

end