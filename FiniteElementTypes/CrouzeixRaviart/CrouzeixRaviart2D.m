classdef CrouzeixRaviart2D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        normalVectors
        n_nodes
        nodes
        barycentricCoords
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = CrouzeixRaviart2D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            syms x y
            hold on
            for i=1:obj.n_vertices
                subplot(1,3,i)
                func = piecewise(x+y<=1,obj.shapeFunctions{i});
                fsurf(func,[0 1])
                xlabel('x')
                ylabel('y')
                zlabel('z')
                title("i:"+i)
                grid on
            end
            hold off
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.computeVertices()
            obj.computeNormalVectors()
            obj.computeBarycentricCoords()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = [sqrt(2)/2,sqrt(2)/2;0,-1;-1,0];
        end
        
        function computeBarycentricCoords(obj)
            syms x y
            for i = 1:obj.n_vertices
                if i~=obj.n_vertices
                    obj.barycentricCoords{i} = symfun(1-dot([x,y]-obj.vertices(i,:),obj.normalVectors(i,:))/dot(obj.vertices(i+1,:)-obj.vertices(i,:),obj.normalVectors(i,:)),[x,y]);
                else
                    obj.barycentricCoords{i} = symfun(1-dot([x,y]-obj.vertices(i,:),obj.normalVectors(i,:))/dot(obj.vertices(1,:)-obj.vertices(i,:),obj.normalVectors(i,:)),[x,y]);
                end
            end
        end
        
        function computeShapeFunctions(obj)
            for i = 1:obj.n_vertices
                obj.shapeFunctions{i} = 1-2*obj.barycentricCoords{i};
            end
        end
        
    end
    
end