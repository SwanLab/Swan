classdef PlotterDensity < handle
    
    properties (Access = private)
        patchHandle
    end

    properties (Access = private)
        density
    end

    methods (Access = public)
        
        function obj = PlotterDensity(cParams)
            obj.createFigure(cParams.mesh);
        end
        
        function plot(obj,rho)
            funp0   = rho.project('P0');
            rhoElem = squeeze(funp0.fValues);
            set(obj.patchHandle,'FaceVertexAlphaData',rhoElem,'FaceAlpha','flat'); 
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
        end
        
        function createFigure(obj,mesh)
            m = mesh;
            figure;
            set(gcf,'Pointer','arrow','NumberTitle','off');
            hold on
            axis off
            axis equal
            axes = gcf().Children;
            obj.patchHandle = patch(axes,'Faces',m.connec,'Vertices',m.coord,...
                'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
        end
            
    end
    
end