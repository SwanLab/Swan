classdef LineMesh < Mesh
    
    properties (Access = public)
        geometryType = 'Line';
    end
    
    properties (Access = private)
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = LineMesh(cParams)
            obj = obj@Mesh(cParams);
            % obj.init(cParams)
        end

        function detJ = computeJacobianDeterminant(obj,xV)
            J = obj.computeJacobian(xV);
            detJ = squeeze(pagenorm(J));
        end
        
        function invJ = computeInverseJacobian(obj,xV)
            detJ = obj.computeJacobianDeterminant(xV);
            invJ = 1./detJ;
        end
        
        function plot(obj)
            p = patch('vertices',obj.coord,'faces',obj.connec);
            p.EdgeColor = 'b';
            p.EdgeAlpha = 1;
            p.EdgeLighting = 'flat';
            p.LineWidth = 0.5;
            p.LineStyle = '-';
            axis('equal');
            nodes = unique(obj.connec(:));
            if size(obj.coord,2) == 3
                x = obj.coord(:,1);
                y = obj.coord(:,2);
                z = obj.coord(:,3);
                hold on
                p = plot3(x(nodes),y(nodes),z(nodes),'.r');
                p.MarkerSize = 6;
            else
                x = obj.coord(:,1);
                y = obj.coord(:,2);
                hold on
                p = plot(x(nodes),y(nodes),'.r');
                p.MarkerSize = 14;
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
        end
        
    end
    
end