classdef LagrangeTensorProduct1D < handle
   
    properties (Access = private)
        k
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        n_nodes
        nodes
        lagrangePolynomials
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = LagrangeTensorProduct1D(k)
            obj.init(k);
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            hold on
            for i = 1:obj.k+1
                subplot(1,obj.k+1,i)
                fplot(obj.shapeFunctions{1,i},[0 1])
                xlabel('x')
                ylabel('y')
                title("i:"+i)
                grid on
            end
            hold off
        end
        
    end
    
    
    methods (Access = private)
       
        function init(obj,k)
            obj.k = k;
            obj.computeVertices()
            obj.computeNodes()
            obj.computeLagrangePolinomyals()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 2;
            obj.vertices = [0,1];
        end
        
        function computeNodes(obj)
            obj.n_nodes = obj.k+1;
            
            if obj.k == 1
                obj.nodes.coord = [0,1];
                obj.nodes.index = [1,2];
            elseif obj.k == 2
                obj.nodes.coord = [0,0.5,1];
                obj.nodes.index = [1,0;1,2;2,0];
            elseif obj.k == 3
                obj.nodes.coord = [0,1/3,2/3,1];
                obj.nodes.index = [1,0;1,2;2,1;2,0];
            end
        end
        
        function computeLagrangePolinomyals(obj)
            syms x
            for i = 1:obj.n_nodes
                func = 1;
                for j = 1:obj.n_nodes
                    if i~=j
                        func = func*(x-obj.nodes.coord(j))/(obj.nodes.coord(i)-obj.nodes.coord(j));
                    end
                end
                obj.lagrangePolynomials{i} = func;
            end
        end
        
        function computeShapeFunctions(obj)
            obj.shapeFunctions = obj.lagrangePolynomials;
        end
        
    end
    
end