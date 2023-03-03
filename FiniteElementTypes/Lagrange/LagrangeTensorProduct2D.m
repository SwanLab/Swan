classdef LagrangeTensorProduct2D < handle
   
    properties (Access = public)
        polynomialOrder
        n_vertices
        vertices
        ndofs
        nodes
        lagrangePolynomials
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = LagrangeTensorProduct2D(polynomialOrder)
            obj.init(polynomialOrder);
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
            
        end
        
    end
    
    
    methods (Access = private)
       
        function init(obj,polynomialOrder)
            obj.polynomialOrder = polynomialOrder;
            obj.computeVertices();
            obj.computeNdof();
            obj.computeNodes();
            obj.computeLagrangePolinomyals();
            obj.computeShapeFunctions();
        end
        
        function computeVertices(obj)
            obj.n_vertices = 4;
            obj.vertices = [0,0;0,1;1,0;1,1];
        end
        
        function computeNdof(obj)
            obj.ndofs = (k+1)^2;
        end
        
        function computeNodes(obj)
            k = obj.polynomialOrder;
            node = zeros(k+1,k+1,2);
            for i = 1:k+1
                for j = 1:k+1
                    node(i,j,1)=(i-1)/k;
                    node(i,j,2)=(j-1)/k;
                end
            end
            obj.nodes = node;
        end
        
        function computeLagrangePolinomyals(obj)
            k = obj.polynomialOrder;
            syms x y
            lagrangePolynomial = cell(2,k+1);
            for i = 1:(k+1)
                func1 = 1;
                func2 = 1;
                for j = 1:(k+1)
                    if i~=j
                        func1 = func1*(x-obj.nodes(j,i,1))/(obj.nodes(i,i,1)-obj.nodes(j,i,1));
                        func2 = func2*(y-obj.nodes(i,j,2))/(obj.nodes(i,i,2)-obj.nodes(i,j,2));
                    end
                end
                lagrangePolynomial{1,i} = func1;
                lagrangePolynomial{2,i} = func2;
            end
            obj.lagrangePolynomials = lagrangePolynomial;
        end
        
        function computeShapeFunctions(obj)
            k = obj.polynomialOrder;
            shapeFunc = cell(k+1);
            for i = 1:(k+1)
                for j = 1:(k+1)
                    shapeFunc{i,j} = matlabFunction(obj.lagrangePolynomials{1,i}*obj.lagrangePolynomials{2,j});
                end 
            end
            obj.shapeFunctions = cell(k+1);
        end
        
    end
    
end