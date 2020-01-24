classdef GradientVariationWithThetaComputer < handle
    
    properties (Access = private)
        levelSet
        radius
        rPerimeter
        mesh
        domainLength
        filePlotName        
        circunferenceMesh
        
        gradientCircunf
        theta
        nEpsilon
        nCell
        legendPlot
    end
    
    methods (Access = public)
        
        function obj = GradientVariationWithThetaComputer(cParams)
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
            obj.circunferenceMesh = obj.computeCircunferenceMesh();
            obj.nEpsilon        = size(obj.rPerimeter.epsilons,2);
            obj.nCell           = size(obj.circunferenceMesh.connec,1);
            obj.theta           = zeros(obj.nCell,obj.nEpsilon);           
            obj.gradientCircunf = zeros(obj.nCell,obj.nEpsilon);                        
        end
        
        function cMesh = computeCircunferenceMesh(obj)
            s = obj.createMeshUnfittedSettings();
            cMesh = UnfittedMesh(s);
            cMesh.compute(obj.levelSet.value);
        end        
        
        function s = createMeshUnfittedSettings(obj)
            mBackground   = obj.levelSet.mesh;
            interpolation = Interpolation.create(mBackground,'LINEAR');
            sM.unfittedType            = 'BOUNDARY';
            sM.meshBackground          = mBackground;
            sM.interpolationBackground = interpolation;
            sM.includeBoxContour       = false;
            s = SettingsMeshUnfitted(sM);
        end
        
        function computeGradientInCircunference(obj,iepsilon)
            gradient = obj.rPerimeter.perimetersGradient(:,iepsilon);
            for icell = 1:obj.nCell
                gNodalInCell  = obj.computeFnodalInCell(icell,gradient);
                xPos          = obj.computeXpos(icell);
                gInXpos       = obj.interpolateValue(xPos,gNodalInCell);
                obj.gradientCircunf(icell,iepsilon) = gInXpos;
            end
        end        
        
        function computeTheta(obj,iepsilon)
            cMesh = obj.circunferenceMesh;
            for icell = 1:obj.nCell
                nodeA  = cMesh.connec(icell,1);
                nodeB  = cMesh.connec(icell,2);
                coordA = cMesh.coord(nodeA,:);
                coordB = cMesh.coord(nodeB,:);
                yA = coordA(:,2)-0.5;
                xA = coordA(:,1)-0.5;
                thet = atan2(yA,xA)*180/pi;
                obj.theta(icell,iepsilon) = thet;
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
            p{1} = plot([-180 180],[(-1/r) (-1/r)]);
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
            [n,d] = rat(h/L);
            tit = ['$h/L = ',num2str(n),'/',num2str(d),'$'];
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
        
        function fNodalInCell = computeFnodalInCell(obj,icell,fNodal)
            cellGlobal   = obj.circunferenceMesh.cellContainingSubcell(icell);
            nodeTriangle = obj.mesh.connec(cellGlobal,:);
            fNodalInCell = fNodal(nodeTriangle(:),:);
        end
        
        function x1pos = computeXpos(obj,icell)
            x1pos = obj.circunferenceMesh.subcellIsoCoords(icell,1,:);
            x1pos = squeeze(x1pos);
        end
        
        function fInterp = interpolateValue(obj,xpos,fNodal)
            interpolation = Interpolation.create(obj.mesh,'LINEAR');
            interpolation.computeShapeDeriv(xpos);
            shape = interpolation.shape;
            fInterp = 0;
            for inode = 1:length(shape)
                fInterp = fInterp + fNodal(inode)*shape(inode);
            end
        end
        
    end
    
  
end