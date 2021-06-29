classdef PlotterDensity < handle
    
    properties (Access = private)
       filter 
       patchHandle
       axes
       figHandle
       nFigures
    end
    
    properties (Access = private)
       figureNumber
       designVariable
    end
    
    methods (Access = public)
        
        function obj = PlotterDensity(cParams)
            obj.init(cParams)
        end
        
        function plot(obj,cost,t,beta,incX)         
            obj.plotDensity();
            obj.plotCost(cost);
            obj.plotLineSearch(t);   
            obj.plotBeta(beta);
            obj.plotIncX(incX);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)            
            obj.designVariable = cParams.designVariable;
            obj.nFigures = 5;
            obj.createFigure();
            obj.createFilter();
            obj.createAxisAndPatchHandle();
        end    
        
        function createFigure(obj)
            fh = figure('units', 'pixels');        
            fh.set('Position',[4000 1500 3000 500])
            obj.figureNumber = fh;           
        end        
    
        function createAxisAndPatchHandle(obj)
            figure(obj.figureNumber)
            subplot(1,obj.nFigures,5)
            obj.figHandle = obj.figureNumber;
            set(obj.figHandle,'Pointer','arrow','NumberTitle','off');
            title('Density')
            hold on
            axis off
            axis equal
            obj.axes = obj.figHandle.Children;%CurrentAxes;
            mesh = obj.designVariable.mesh;
            obj.patchHandle = patch(obj.axes,'Faces',mesh.connec,'Vertices',mesh.coord,...
                'FaceAlpha','flat','EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
            set(obj.axes,'ALim',[0, 1],'XTick',[],'YTick',[]);
            obj.figHandle.Color = 'white';
            colormap([0 0 0]);
        end
        
        function createFilter(obj)
            s.filterType = 'P1';
            s.domainType = 'INTERIOR';
            s.designVar = obj.designVariable;
            s.quadratureOrder = 'LINEAR';
            s = SettingsFilter(s);
            s.femSettings.scale = 'MACRO';
            filterP1 = Filter_P1_Density(s);
            filterP1.preProcess();
            obj.filter = filterP1;
        end
        
        function plotDensity(obj)
            rho = obj.designVariable.value;
            rhoElem = obj.filter.getP0fromP1(rho);
            set(obj.patchHandle,'FaceVertexAlphaData',rhoElem,'FaceAlpha','flat'); 
            drawnow            
        end
        
        function plotCost(obj,cost)
            figure(obj.figureNumber)
            subplot(1,obj.nFigures,1)
            plot(cost);
            title('Cost')            
        end
        
        function plotLineSearch(obj,t)
            figure(obj.figureNumber)            
            subplot(1,obj.nFigures,2)
            plot(t);
            title('LineSearch')              
        end
        
        function plotBeta(obj,beta)
            figure(obj.figureNumber)            
            subplot(1,obj.nFigures,3)
            plot(beta);
            title('$\beta$','Interpreter','Latex');                          
        end
        
        function plotIncX(obj,incX)
            figure(obj.figureNumber)            
            subplot(1,obj.nFigures,4)
            plot(incX);
            title('$\Delta x$','Interpreter','Latex')                                                  
        end        
    
    end
    
end
