classdef ExteriorGradientVariationWithThetaComputer < handle
    
    properties (Access = private)
        levelSet
        radius
        rPerimeter
        mesh
        domainLength
        filePlotName  
        exteriorNodes
        
        gradientCircunf
        theta
        nEpsilon
        nNode
        nodes
        legendPlot        
    end
    
    methods (Access = public)
        
        function obj = ExteriorGradientVariationWithThetaComputer(cParams)
            obj.init(cParams)
        end
            
        function compute(obj)            
            for iepsilon = 1:obj.nEpsilon
                obj.computeGradientInCircunference(iepsilon);
                obj.computeTheta(iepsilon);
                obj.computeLegend(iepsilon);
            end
            obj.plotGradientVsTheta();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filePlotName = cParams.filePlotName;
            obj.mesh         = cParams.mesh;
            obj.levelSet     = cParams.levelSet;
            obj.radius       = cParams.radius;
            obj.rPerimeter   = cParams.regularizedPerimeter;
            obj.domainLength = cParams.domainLength;
            obj.exteriorNodes = obj.computeExteriorNodes();
            obj.nEpsilon        = size(obj.rPerimeter.epsilons,2);
            obj.nNode           = size(obj.exteriorNodes,1);
            obj.theta           = zeros(obj.nNode,obj.nEpsilon);           
            obj.gradientCircunf = zeros(obj.nNode,obj.nEpsilon);                        
        end
        
        function eNodes = computeExteriorNodes(obj)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            isExterior = abs(sqrt(x.^2 + y.^2) - obj.radius) < 1e-6;
            allNodes = 1:length(obj.mesh.coord(:,1));
            eNodes(:,1) = allNodes(isExterior);
        end
        
        function computeGradientInCircunference(obj,iepsilon)
            gradient = obj.rPerimeter.perimetersGradient(:,iepsilon);
            for inode = 1:obj.nNode
                node            = obj.exteriorNodes(inode);
                gInExteriorNode = gradient(node,1);
                obj.gradientCircunf(inode,iepsilon) = gInExteriorNode;
            end
        end        
        
        function computeTheta(obj,iepsilon)
            for inode = 1:obj.nNode
                node = obj.exteriorNodes(inode);
                coord = obj.mesh.coord(node,:);
                yA = coord(:,2);
                xA = coord(:,1);
                thet = atan2(yA,xA)*180/pi;
                obj.theta(inode,iepsilon) = thet;
            end
        end
        
        function computeLegend(obj,iepsilon)
            epsilons = obj.rPerimeter.epsilons;
            h = obj.mesh.computeMeanCellSize;
            epsStr = num2str(epsilons(iepsilon)/h);
            legStr = ['$dPer^R_\varepsilon \ \textrm{with} \ \varepsilon/h \, = \,',epsStr,'$'];
            obj.legendPlot{end+1} = legStr;
        end
        
        function plotGradientVsTheta(obj)
            f = figure();
            r = obj.radius;
            p{1} = plot([-180 180],[(1/r) (1/r)]);
            leg{1} = '$\kappa = 1/R$';
            hold on
            for iepsilon = 2:obj.nEpsilon
                [x,y] = obj.obtainDataByIncreasingTheta(iepsilon);
                p{end+1} = plot(x,y,'+-');
                leg{end+1} = obj.legendPlot{iepsilon};
            end
            legObj = legend(leg);
            set(legObj,'Interpreter','latex','Location','Best');
            xlabel('$\theta$','interpreter','latex')
            set(gca, 'XTickLabel',[])      
            xticks([-pi -pi/2 0 pi/2 pi]*180/pi) 
            xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
            title(obj.computeTitle,'interpreter','latex')
            obj.printPlot(f,p);
        end
        
        function tit = computeTitle(obj)
            h = obj.mesh.computeMeanCellSize;
            L = obj.domainLength;
            %tit = ['$h/L = ',num2str(n),'/',num2str(d),'$'];
            tit = ['$h/L = ',num2str(round(h/L,3)),'$'];                                    
        end
        
        function [x,y] = obtainDataByIncreasingTheta(obj,iepsilon)
            [~,isort] = sort(obj.theta(:,iepsilon));
            x = obj.theta(isort,iepsilon);
            y = obj.gradientCircunf(isort,iepsilon); 
        end    
        
        function printPlot(obj,f,h)
            outputName = [obj.filePlotName,'GradientVsTheta'];
            printer = plotPrinter(f,h);
            printer.print(outputName);
        end                
        
    end
    
  
end