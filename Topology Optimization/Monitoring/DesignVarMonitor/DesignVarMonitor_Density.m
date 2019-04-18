classdef DesignVarMonitor_Density < DesignVarMonitor_Abstract
    
    properties (Access = protected)
        designVarName = 'Density - \rho';
    end
    
    properties (Access = private)
        rhoElem
        filter
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_Density(designVar,showBC)
            obj@DesignVarMonitor_Abstract(designVar,showBC);
            obj.createFilter(designVar);
        end
        
        function plot(obj)
            obj.filterDensity();
            set(obj.patchHandle,'FaceVertexAlphaData',obj.rhoElem);
        end
        
    end
    
    methods (Access = protected)
        
        function initPlotting(obj)
            obj.patchHandle = patch(obj.axes,'Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,...
                'FaceAlpha','flat','EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
            set(gca,'ALim',[0, 1],'XTick',[],'YTick',[]);
            
            obj.BCplotter.plot();
        end
        
    end
    
    methods (Access = protected, Static)
        
        function color = getColor()
            color = [0 0 0];
        end
        
    end
    
    methods (Access = private)
        
        function createFilter(obj,designVar)
            filterParams = SettingsFilter();
            filterParams.filterType = 'P1';
            filterParams.domainType = 'INTERIOR';
            filterParams.designVar = designVar;
            obj.filter = Filter_P1_Density(filterParams);
            obj.filter.preProcess();
        end
        
        function rhoElem = filterDensity(obj)
            rho = obj.designVar.value;
            if obj.isNodal(rho)
                rhoElem = obj.filter.getP0fromP1(rho);
            elseif obj.isElemental(rho)
                rhoElem = rho;
            else
                error('Invalid density vector size')
            end
            obj.rhoElem = double(rhoElem);
        end
        
        function itIs = isNodal(obj,x)
            n = size(x,1);
            nnode = size(obj.mesh.coord,1);
            itIs = n == nnode;
        end
        
        function itIs = isElemental(obj,x)
            n = size(x,1);
            nelem = size(obj.mesh.connec,1);
            itIs = n == nelem;
        end
        
        
    end
    
end