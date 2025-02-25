classdef PlotterDensity < handle
    
    properties (Access = private)
        patchHandle
    end

    properties (Access = private)
        mesh
        designVariable
    end

    methods (Access = public)
        
        function obj = PlotterDensity(cParams)
            obj.init(cParams);
            obj.createFigure();
        end
        
        function plot(obj)
            rho     =   obj.designVariable.fun;
            funp0   = rho.project('P0');
            rhoElem = squeeze(funp0.fValues);
            set(obj.patchHandle,'FaceVertexAlphaData',rhoElem,'FaceAlpha','flat'); 
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.designVariable = cParams.designVariable;
        end
        
        function createFigure(obj)
            figure;
            set(gcf,'Pointer','arrow','NumberTitle','off');
            hold on
            axis off
            axis equal
            axes = gcf().Children;
            obj.patchHandle = patch(axes,'Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,...
                'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
        end
            
    end
    
end