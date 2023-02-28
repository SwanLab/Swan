classdef LagrangeSimplicial2D < handle
   
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
    
        function obj = LagrangeSimplicial2D(k)
            obj.init(k);
        end
        
        function plotShapeFunctions(obj)
            s.coord = obj.vertices;
            s.connec = [1 2 3];
            m = Mesh(s);
            
            obj.fig = figure();
            for i=1:obj.n_nodes
                subplot(obj.k+1,obj.k+1,i)
                m.plot()
                
                points = obj.evaluateFunction(m,i);
                plot3(points(:,1),points(:,2),points(:,3),'.');
                
                zlim([-0.5 1])
                xlabel('x'); ylabel('y'); zlabel('z')
                title("i:"+string(i-1))
                grid on
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj,k)
            obj.k = k;
            obj.computeVertices()
            obj.computeNodes()
            obj.computeShapeFunctions()
            obj.storeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeNodes(obj)
            obj.n_nodes = nchoosek(2+obj.k,obj.k);
            
            obj.nodes = zeros(obj.k+1,obj.k+1,2);
            for i = 1:obj.k+1
                for j = 1:obj.k+1
                    if (i+j)<=(obj.k+2)
                        obj.nodes(i,j,1)=(i-1)/obj.k;
                        obj.nodes(i,j,2)=(j-1)/obj.k;
                    end
                end
            end
        end
        
        function computeShapeFunctions(obj)
            syms x y
            X = obj.assemblyBasis();
            
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    if (i+j)<=(obj.k+2)
                        s = obj.computeShapeFunctionCoefficients(X,i,j);
                        I = obj.computeSimplexIndeces(i,j);
                        obj.shapeFunctions{I} = X*s;
                    end
                end
            end
        end
        
        function s = computeShapeFunctionCoefficients(obj,X,i,j)
            A = obj.assemblyShapeFunctionCoefficientsLHS(X);
            B = obj.assemblyShapeFunctionCoefficientsRHS(i,j);
            s = A\B;
        end
        
        function storeShapeFunctions(obj)
            syms x y
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    if (i+j)<=(obj.k+2)
                        I = obj.computeSimplexIndeces(i,j);
                        f = matlabFunction(obj.shapeFunctions{I},'Vars',[x y]);
                        obj.shapeFunctions{I} = f;
                    end
                end
            end
        end
        
        function X = assemblyBasis(obj)
            syms x y
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    if (i+j)<=(obj.k+2)
                        I = obj.computeSimplexIndeces(i,j);
                        X(I) = x^(i-1)*y^(j-1);
                    end
                end
            end
        end
        
        function B = assemblyShapeFunctionCoefficientsRHS(obj,i,j)
            I = obj.computeSimplexIndeces(i,j);
            B = zeros(obj.n_nodes,1);
            B(I) = 1;
        end
        
        function A = assemblyShapeFunctionCoefficientsLHS(obj,X)
            syms x y
            A = zeros(obj.n_nodes);
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    if (i+j)<=(obj.k+2)
                        I = obj.computeSimplexIndeces(i,j);
                        node(:) = obj.nodes(i,j,:);
                        A(I,:) = subs(X,[x y],node);
                    end
                end
            end
        end
        
        function I = computeSimplexIndeces(obj,i,j)
            n = -(obj.k+2);
            for m = 1:i
                n = n + obj.k + 1 - (m-2);
            end
            I = n+j;
        end
        
        function inside_points = evaluateFunction(obj,m,I)
            A = m.coord(1,:);
            B = m.coord(2,:);
            C = m.coord(3,:);

            points = rand(10000,2);

            points(:,1) = points(:,1) * max([B(1) C(1)]) + min([A(1) B(1) C(1)]);
            points(:,2) = points(:,2) * max([C(2) sqrt(3)/2]) + min([A(2) B(2) C(2)]);

            inside_points = [];
            for i = 1:size(points,1)
                P = points(i,:);
                if P(1)+P(2) <= 1
                    inside_points = [inside_points; P];
                end
            end
            inside_points(:,3) = obj.shapeFunctions{I}(inside_points(:,1),inside_points(:,2));
        end
        
    end
    
end