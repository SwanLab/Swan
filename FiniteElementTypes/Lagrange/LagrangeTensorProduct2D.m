classdef LagrangeTensorProduct2D < handle
   
    properties (Access = public)
        k
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
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            set(groot,'defaultAxesTickLabelInterpreter','latex');
            obj.fig = figure();
            hold on
            for i = 1:obj.k+1
                for j = 1:obj.k+1
                    subplot(obj.k+1,obj.k+1,(obj.k+1)*(i-1)+j)
                    fsurf(obj.shapeFunctions{i,j},[0 1 0 1])
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    title("i:"+string(i-1)+", j:"+string(j-1))
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
            
            obj.nodes = zeros(obj.k+1,obj.k+1,2);
            for i = 1:obj.k+1
                for j = 1:obj.k+1
                    obj.nodes(i,j,1)=(i-1)/obj.k;
                    obj.nodes(i,j,2)=(j-1)/obj.k;
                end
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
                        func1 = func1*(x-obj.nodes(j,i,1))/(obj.nodes(i,i,1)-obj.nodes(j,i,1));
                        func2 = func2*(y-obj.nodes(i,j,2))/(obj.nodes(i,i,2)-obj.nodes(i,j,2));
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
                    obj.shapeFunctions{i,j} = matlabFunction(obj.lagrangePolynomials{1,i}*obj.lagrangePolynomials{2,j});
                end 
            end
        end
        
    end
    
end