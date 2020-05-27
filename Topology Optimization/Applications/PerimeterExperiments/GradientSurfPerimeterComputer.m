classdef GradientSurfPerimeterComputer < handle
    
    properties (Access = private)
        outPutFolder
        mesh
        rPerimeter
        iMesh
        domainLength
        figureID
    end
    
    methods (Access = public)
        
        function obj = GradientSurfPerimeterComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            epsilons = obj.rPerimeter.epsilons;
            nEpsilon = length(epsilons);
            for iEpsilon = 1:nEpsilon
                obj.plotSurf(iEpsilon);
                obj.addXYlabel();
                obj.addTitle(iEpsilon);
                obj.printSurf(iEpsilon);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.outPutFolder = cParams.outPutFolder;
            obj.mesh         = cParams.mesh;
            obj.rPerimeter   = cParams.rPerimeter;
            obj.iMesh        = cParams.iMesh;
            obj.domainLength = cParams.domainLength;
        end
        
        function addXYlabel(obj)
            xlabel('$x$','interpreter','latex')
            ylabel('$y$','interpreter','latex')
        end
        
        function plotSurf(obj,iEpsilon)
            coord = obj.mesh.coord;
            x = coord(:,1);
            y = coord(:,2);
            z = obj.rPerimeter.perimetersGradient(:,iEpsilon);
            tri = delaunay(x,y);
            f = figure();
            trisurf(tri,x,y,z);
            shading interp
            obj.figureID = f;
        end
        
        function printSurf(obj,iEpsilon)
            part1 = [obj.outPutFolder,'GradientMesh',num2str(obj.iMesh)];
            part2 = ['Epsilon',num2str(iEpsilon)];
            outFile = [part1,part2];
            printer = surfPrinter(obj.figureID);
            printer.print(outFile)
        end
        
        function addTitle(obj,iEpsilon)
            epsilons = obj.rPerimeter.epsilons;
            h = obj.mesh.computeMeanCellSize;
            epsStr = num2str(epsilons(iEpsilon)/h);
            firstStr = 'dPer^R_\varepsilon \ \textrm{with} \ \,';
            secondStr = '\varepsilon/h \, = \, ';
            thirdStr = '\ \textrm{and} \ h/L \, = \, ';
            L = obj.domainLength;
            [n,d] = rat(h/L);
            hStr = [num2str(n),'/',num2str(d)];
            titleStr = ['$',firstStr,secondStr,epsStr,thirdStr,hStr,'$'];
            title(titleStr,'interpreter','latex');
        end        
        
    end
    
end