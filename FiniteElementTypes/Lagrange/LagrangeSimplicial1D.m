classdef LagrangeSimplicial1D < handle
   
    properties (Access = public)
        k
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
    
        function obj = LagrangeSimplicial1D(k)
            obj.init(k);
        end
    
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            set(groot,'defaultAxesTickLabelInterpreter','latex');
            obj.fig = figure();
            hold on
            for i = 1:obj.k+1
                subplot(1,obj.k+1,i)
                fplot(obj.shapeFunctions{i},[0 1])
                xlabel('x')
                ylabel('y')
                title("$i="+string(i-1)+"$")
                grid on
            end
            hold off
        end
        
    end
    
    
    methods (Access = private)
       
        function init(obj,k)
            obj.k = k;
            obj.computeVertices()
            obj.computeNormalVectors()
            obj.computeNodes()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 2;
            obj.vertices = [0,1];
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = [1,0];
        end
        
        function computeNodes(obj)
            obj.n_nodes = nchoosek(1+obj.k,obj.k);
            for i=1:obj.n_nodes
                obj.nodes.coord(i) = (i-1)/obj.k;
                obj.nodes.index(i) = i;
            end
        end

        function computeShapeFunctions(obj)
            for i = 1:obj.n_nodes
                syms x
                shapeFunc = 1;
                for j = 1:obj.n_nodes
                    if i~=j
                        shapeFunc = shapeFunc * (x-obj.nodes.coord(j))/(obj.nodes.coord(i)-obj.nodes.coord(j));
                    end
                end
                obj.shapeFunctions{i} = matlabFunction(shapeFunc);
            end
        end
        
    end
    
end