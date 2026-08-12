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
            caxis([-0.8 0.8])
        end
        % function plot(obj)
        %     rho = obj.designVariable.fun;
        %     funp0 = rho.project('P0');
        %     rhoElem = squeeze(funp0.fValues);
        %     set(obj.patchHandle, 'FaceVertexCData', rhoElem, 'FaceColor', 'flat');
        %     caxis([-0.6 0.6])
        %     colormap(gray);  
        %     colorbar;
        % end
        
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
        % function createFigure(obj)
        %     figure;
        %     set(gcf,'Pointer','arrow','NumberTitle','off');
        %     hold on
        %     axis off
        %     axis equal
        %     axes = gcf().Children;
        %     obj.patchHandle = patch(axes,'Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,...
        %         'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
        %     colormap(gray);  % escala de cinza (preto e branco)
        %     colorbar
        % end

    end
    methods (Access = private, Static)

        function c = redtoblue(m)
            if nargin < 1; m = 64; end
            half = floor(m/2);
            r1 = linspace(0, 1, half)';
            g1 = linspace(0, 0, half)';
            b1 = linspace(1, 0, half)';
            r2 = linspace(1, 0, m-half)';
            g2 = linspace(0, 0, m-half)';
            b2 = linspace(0, 1, m-half)';
            c  = [r1, g1, b1; r2, g2, b2];
        end

    end
    
end