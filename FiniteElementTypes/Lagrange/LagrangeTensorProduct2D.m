classdef LagrangeTensorProduct2D < handle
   
    properties (Access = private)
        k
    end

    
    properties (Access = private)
        n_vertices
        vertices
        n_nodes
        nodes
        lagrangePolynomials
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = LagrangeTensorProduct2D(k)
            obj.init(k);
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            hold on
            for i = 1:obj.k+1
                for j = 1:obj.k+1
                    subplot(obj.k+1,obj.k+1,(obj.k+1)*(i-1)+j)
                    fsurf(obj.shapeFunctions{i,j},[0 1 0 1])
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    title("i:"+i+", j:"+j)
                end
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
            obj.n_vertices = 4;
            obj.vertices = [0,0;0,1;1,0;1,1];
        end
        
        function computeNodes(obj)
            obj.n_nodes = (obj.k+1)^2;
            
            if obj.k == 1
                obj.nodes.coord = [0,0;0,1;1,0;1,1];
                obj.nodes.index = [1,1;1,2;2,1;2,2];
            elseif obj.k == 2
                obj.nodes.coord = [0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1];
                obj.nodes.index = [1,1;1,2;1,3;2,1;2,2;2,3;3,1;3,2;3,3];
            elseif obj.k == 3
                obj.nodes.coord = [0,0;0,1/3;0,2/3;0,1;1/3,0;1/3,1/3;1/3,2/3;1/3,1;2/3,0;2/3,1/3;2/3,2/3;2/3,1;1,0;1,1/3;1,2/3;1,1];
                obj.nodes.index = [1,1;1,2;1,3;1,4;2,1;2,2;2,3;2,4;3,1;3,2;3,3;3,4;4,1;4,2;4,3;4,4];
            end
        end
        
        function computeLagrangePolinomyals(obj)
            syms x y
            obj.lagrangePolynomials = cell(2,obj.k+1);
            for i = 1:(obj.k+1)
                func1 = 1;
                func2 = 1;
                for j = 1:(obj.k+1)
                    if i~=j
                        func1 = func1*(x-obj.nodes.coord((j-1)*(obj.k+1)+1,1))/(obj.nodes.coord((i-1)*(obj.k+1)+1,1)-obj.nodes.coord((j-1)*(obj.k+1)+1,1));
                        func2 = func2*(y-obj.nodes.coord(j,2))/(obj.nodes.coord(i,2)-obj.nodes.coord(j,2));
                    end
                end
                obj.lagrangePolynomials{1,i} = func1;
                obj.lagrangePolynomials{2,i} = func2;
            end
        end
        
        function computeShapeFunctions(obj)
            obj.shapeFunctions = cell(obj.k+1);
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    obj.shapeFunctions{i,j} = obj.lagrangePolynomials{1,i}*obj.lagrangePolynomials{2,j};
                end 
            end
        end
        
    end
    
end