classdef GradientSurfPerimeterComputer < handle

    properties (Access = private)
        outPutFolder
        iMesh
        figureID
        gExperiment
    end

    properties (Access = private)
        inputFiles
        outputFolder
        levelSetParams
        circleCase
    end

    methods (Access = public)

        function obj = GradientSurfPerimeterComputer(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            for im = 1:numel(obj.inputFiles)
                obj.iMesh = im;
                obj.createGradientVariationExperiment();
                nEpsilon = length(obj.gExperiment.regularizedPerimeter.epsilons);
                for iEpsilon = 1:nEpsilon
                    obj.plotSurf(iEpsilon);
                    obj.addXYlabel();
                    obj.addTitle(iEpsilon);
                    obj.printSurf(iEpsilon);
                end
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.inputFiles   = cParams.inputFiles;
            obj.outputFolder = cParams.outputFolder;
            obj.levelSetParams = cParams.levelSetParams;
            obj.circleCase = cParams.circleCase;
        end

        function createGradientVariationExperiment(obj)
            s.levelSetParams = obj.levelSetParams;
            s.inputFile      = obj.inputFiles{obj.iMesh};
            s.iMesh          = obj.iMesh;
            g = GradientVariationExperiment(s);
            obj.gExperiment = g;
        end

        function addXYlabel(obj)
            xlabel('$x$','interpreter','latex')
            ylabel('$y$','interpreter','latex')
        end

        function plotSurf(obj,iEpsilon)
            coord = obj.gExperiment.backgroundMesh.coord;
            x = coord(:,1);
            y = coord(:,2);
            dPer = obj.gExperiment.regularizedPerimeter.perimetersGradient;
            z = dPer(:,iEpsilon);
            tri = delaunay(x,y);
            f = figure();
            h = trisurf(tri,x,y,z);
            shading interp
            obj.figureID = f;
            ar = get(gca,'DataAspectRatio');
            obj.plotBoundaryCutMesh(tri,x,y,z)
            set(gca,'DataAspectRatio',ar);
            view(-60,25);
        end

        function plotBoundaryCutMesh(obj,tri,x,y,z)
            m = obj.computeBoundaryCutMesh();
            %nodes = unique(m.connec(:));
            Xi = m.coord(:,1);
            Yi = m.coord(:,2);
            Zi = interptri(tri,x,y,z,Xi,Yi);
            s.coord = [Xi Yi Zi];
            %s.connec(:,1) = 1:length(Xi);
            %s.connec(:,2) = 2:length(Xi)+1;
            %s.connec(end,2) = 1;
            %s.coord = m.coord;
            s.connec = m.connec;
            s.kFace = -2;
            mB = Mesh.create(s);
            hold on
            mB.plot
            % axes(axes_h)
            % axis(axis_h)
            % view(a,b)
        end

        function mBCut = computeBoundaryCutMesh(obj)
            g = obj.gExperiment;
            s.backgroundMesh = g.backgroundMesh;
            s.boundaryMesh   = g.boundaryMesh;
            uMesh = UnfittedMesh(s);
            uMesh.compute(g.levelSet);
            bCut = uMesh.boundaryCutMesh;
            mBCut = bCut.mesh;
        end

        function printSurf(obj,iEpsilon)
            part1 = [obj.outputFolder,obj.circleCase,'GradientMesh',num2str(obj.iMesh)];
            part2 = ['Epsilon',num2str(iEpsilon)];
            outFile = [part1,part2];
            printer = surfPrinter(obj.figureID);
            printer.print(outFile)
        end

        function addTitle(obj,iEpsilon)
            epsilons = obj.gExperiment.regularizedPerimeter.epsilons;
            h = obj.gExperiment.backgroundMesh.computeMeanCellSize;
            epsStr = num2str(epsilons(iEpsilon)/h);
            firstStr = 'dPer^R_\varepsilon \ \textrm{with} \ \,';
            secondStr = '\varepsilon/h \, = \, ';
            thirdStr = '\ \textrm{and} \ h/L \, = \, ';
            L = obj.gExperiment.domainLength;
            [n,d] = rat(h/L);
            if n == 1
                hStr = [num2str(n),'/',num2str(d)];
            else
                hStr = num2str(round(h/L,3));
            end
            titleStr = ['$',firstStr,secondStr,epsStr,thirdStr,hStr,'$'];
            title(titleStr,'interpreter','latex');
        end

    end

end