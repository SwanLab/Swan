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
                subplot(obj.k+1,obj.k+1,I)
                m.plot()
                
                points = evaluateFunction(m,I);
                surf(points(:,1),points(:,2),points(:,3));
                
                zlim([-0.5 1])
                xlabel('x'); ylabel('y'); zlabel('z')
                title("I:"+string(I-1))
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
                        I = computeSimplexIndeces(i,j);
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
                        obj.shapeFunctions{i,j} = matlabFunction(obj.shapeFunctions{i,j});
                    end
                end
            end
        end
        
        function X = assemblyBasis(obj)
            syms x y
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    if (i+j)<=(obj.k+2)
                        I = computeSimplexIndeces(i,j);
                        X(I) = x^(i-1)*y^(j-1);
                    end
                end
            end
        end
        
        function B = assemblyShapeFunctionCoefficientsRHS(obj,i,j)
            I = computeSimplexIndeces(i,j);
            B = zeros(obj.n_nodes,1);
            B(I) = 1;
        end
        
        function A = assemblyShapeFunctionCoefficientsLHS(obj,X)
            syms x y
            A = zeros(obj.n_nodes);
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    if (i+j)<=(obj.k+2)
                        I = computeSimplexIndeces(i,j);
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
        
        function points = obj.evaluateFunction(obj,m,I)
            p1 = m.coord(1,:);
            p2 = m.coord(2,:);
            p3 = m.coord(3,:);

            b = rand(1000,2);
            b(:,3) = 1 - b(:,1) - b(:,2);
            points = b(:,1) .* p1 + b(:,2) .* p2 + b(:,3) .* p3;
            points(:,3) = obj.shapeFunctions{I}(points(:,1),points(:,2));
        end
        
    end
    
end