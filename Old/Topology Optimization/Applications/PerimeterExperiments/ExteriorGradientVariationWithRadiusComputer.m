classdef ExteriorGradientVariationWithRadiusComputer < handle

    properties (Access = private)
        mesh
        nameCase
        inputFile
        outputFolder
        rPerimeter
        xNodesToPlot
        gradientInNodesToPlot     
        domainLength
    end
    
    methods (Access = public)
        
        function obj = ExteriorGradientVariationWithRadiusComputer(cParams)
            obj.init(cParams)            
        end
        
        function compute(obj)            
            obj.xNodesToPlot          = obj.computeXinNodesToPlot();
            obj.gradientInNodesToPlot = obj.computeGradientInNodesToPlot();
            obj.plotGradientWithRadius();
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.rPerimeter   = cParams.regularizedPerimeter;
            obj.outputFolder = cParams.outputFolder;
            obj.inputFile    = cParams.inputFile;
            obj.nameCase     = cParams.nameCase;
            obj.domainLength = cParams.domainLength;            
        end

        
        function xN = computeXinNodesToPlot(obj)
            x = obj.mesh.coord(:,1);
            nodesToPlot = obj.obtainNodesToPlot();            
            xN = x(nodesToPlot,1);            
        end
        
        function g = computeGradientInNodesToPlot(obj)
            gradient = obj.rPerimeter.perimetersGradient;
            nodeToPlot = obj.obtainNodesToPlot();
            g = gradient(nodeToPlot,:);
        end
        
        function nodes = obtainNodesToPlot(obj)
            x  = obj.mesh.coord(:,1);
            y  = obj.mesh.coord(:,2);
            x0 = 0;
            y0 = 0; 
            isYequalToY0    = abs(y-y0) < 0.02;
            isXlargerThanX0 = (x0 - x)  < 1e-12;
            allNodes(:,1) = 1:size(obj.mesh.coord,1);
            nodes = allNodes(isYequalToY0 & isXlargerThanX0,1);
            [~,ind] = sort(x(nodes));
            nodes = nodes(ind);
        end        
        
       function plotGradientWithRadius(obj)
            x = obj.xNodesToPlot;
            y = obj.gradientInNodesToPlot;
            epsilons = obj.rPerimeter.epsilons;
            nEpsilon = length(epsilons);
            h = obj.mesh.computeMeanCellSize;
            leg = cell(nEpsilon,1);
            p   = cell(nEpsilon,1);
            f = figure();
            hold on
            for iepsilon = 1:nEpsilon
                p{iepsilon} = plot(x,y(:,iepsilon),'+-');
                epsStr   = num2str(epsilons(iepsilon)/h);
                leg{iepsilon} = ['$dPer^R_\varepsilon \ \textrm{with} \ \varepsilon/h \, = \,',epsStr,'$'];
            end
            legObj = legend(leg);
            set(legObj,'Interpreter','latex','Location','Best');
            xlabel('$r$','interpreter','latex')
            title(obj.computeTitle,'interpreter','latex')            
            obj.printPlot(f,p);        
       end
       
        function tit = computeTitle(obj)
            h = obj.mesh.computeMeanCellSize;
            L = obj.domainLength;
            %[n,d] = rat(h/L);
            %tit = ['$h/L = ',num2str(n),'/',num2str(d),'$'];
            tit = ['$h/L = ',num2str(round(h/L,3)),'$'];                                                
        end       
       
       function printPlot(obj,f,h)
            fileName = strcat(obj.nameCase,obj.inputFile);
            outputName = fullfile(obj.outputFolder,fileName);
            printer = plotPrinter(f,h);
            printer.print(outputName);
        end              
       
    end
    
end