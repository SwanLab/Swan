classdef LagrangeTensorProduct2D < handle
   
    properties (Access = public)
        polynomialOrder
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
            
            k = obj.polynomialOrder;
            s.coord = obj.vertices;
            s.connec = [1 2 3;2 3 4];
            m = Mesh(s);
            for i=1:3
                m = m.remesh(2);
            end
            
%             obj.fig = figure();
            for i = 1:k+1
                for j = 1:k+1
%                     subplot(obj.k+1,obj.k+1,(obj.k+1)*(i-1)+j)
                    figure()
                    m.plot();
                    trisurf(m.connec,m.coord(:,1),m.coord(:,2),obj.shapeFunctions{i,j}(m.coord(:,1),m.coord(:,2)));
                    
                    xlim([0 1]); ylim([0 1]); zlim([-0.5 1]);
                    xlabel('x'); ylabel('y'); zlabel('z');
                    title("i:"+string(i-1)+", j:"+string(j-1))
                    grid on
                end
            end
            hold off
            
        end
        
    end
    
    
    methods (Access = private)
       
        function init(obj,polynomialOrder)
            obj.polynomialOrder = polynomialOrder;
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
            k = obj.polynomialOrder;
            obj.n_nodes = (k+1)^2;
            
            obj.nodes = zeros(k+1,k+1,2);
            for i = 1:k+1
                for j = 1:k+1
                    obj.nodes(i,j,1)=(i-1)/k;
                    obj.nodes(i,j,2)=(j-1)/k;
                end
            end
        end
        
        function computeLagrangePolinomyals(obj)
            k = obj.polynomialOrder;
            syms x y
            obj.lagrangePolynomials = cell(2,k+1);
            for i = 1:(k+1)
                func1 = 1;
                func2 = 1;
                for j = 1:(k+1)
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
            for i = 1:(k+1)
                for j = 1:(k+1)
                    obj.shapeFunctions{i,j} = matlabFunction(obj.lagrangePolynomials{1,i}*obj.lagrangePolynomials{2,j});
                end 
            end
        end
        
    end
    
end