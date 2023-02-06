classdef CrouzeixRaviart1D < handle
   
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
    
        function obj = CrouzeixRaviart1D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            syms x
            hold on
            for i=1:obj.n_vertices
                subplot(1,2,i)
                fplot(obj.shapeFunctions{i},[0 1])
                xlabel('x')
                ylabel('y')
                title("i:"+i)
                grid on
            end
            hold off
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.computeVertices()
            obj.computeBarycentricCoords()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 2;
            obj.vertices = [0,1];
        end
        
        function computeBarycentricCoords(obj)
            syms x
            for i = 1:obj.n_vertices
                if i~=obj.n_vertices
                    obj.barycentricCoords{i} = symfun(1-(x-obj.vertices(i))/(obj.vertices(i+1)-obj.vertices(i)),x);
                else
                    obj.barycentricCoords{i} = symfun(1-(x-obj.vertices(i))/(obj.vertices(1)-obj.vertices(i)),x);
                end
            end
        end
        
        function computeShapeFunctions(obj)
            for i = 1:obj.n_vertices
                obj.shapeFunctions{i} = 1-obj.barycentricCoords{i};
            end
        end
        
    end
    
end