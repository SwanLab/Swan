classdef LagrangeSimplicial2D < handle
   
    properties (Access = public)
        polinomialOrder
        n_vertices
        vertices
        normalVectors
        ndofs
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
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            set(groot,'defaultAxesTickLabelInterpreter','latex');
            
            s.coord = obj.vertices;
            s.connec = [1 2 3];
            m = Mesh(s);
            for i=1:3
                m = m.remesh(2);
            end
            
%             obj.fig = figure();
            for i=1:obj.ndofs
%                 subplot(obj.polinomialOrder+1,obj.polinomialOrder+1,i);
                figure()
                m.plot();
                trisurf(m.connec,m.coord(:,1),m.coord(:,2),obj.shapeFunctions{i}(m.coord(:,1),m.coord(:,2)));
                
                xlim([0 1]); ylim([0 1]); zlim([-0.5 1]);
                xlabel('x'); ylabel('y'); zlabel('z');
                title("i:"+string(i-1));
                grid on
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj,k)
            obj.polinomialOrder = k;
            obj.computeVertices();
            obj.computeNdof();
            obj.computeNodes();
            obj.computeShapeFunctions();
            obj.storeShapeFunctions();
        end

        function computeNdof(obj)
            k = obj.polinomialOrder;                                    
            n = nchoosek(2+k,k);            
            obj.ndofs = n;
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeNodes(obj)
            k = obj.polinomialOrder;                                    
            node = zeros(k+1,k+1,2);
            for i = 1:k+1
                for j = 1:k+1
                    if (i+j)<=(k+2)
                        node(i,j,1)=(i-1)/k;
                        node(i,j,2)=(j-1)/k;
                    end
                end
            end
            obj.nodes = node;
        end
        
        function computeShapeFunctions(obj)
            k = obj.polinomialOrder;
            X = obj.computeBasisInMonomialForm();
            for i = 1:(k+1)
                for j = 1:(k+1)
                    if (i+j)<=(k+2)
                        a = obj.computeShapeFunctionCoefficients(X,i,j);
                        s = obj.computeMonomialIndeces(i,j);
                        obj.shapeFunctions{s} = X*a;
                    end
                end
            end
        end
        
        function s = computeShapeFunctionCoefficients(obj,X,i,j)
            A = obj.assemblyShapeFunctionCoefficientsLHS(X);
            b = obj.assemblyShapeFunctionCoefficientsRHS(i,j);
            s = A\b;
        end
        
        function storeShapeFunctions(obj)
            syms x y
            for i = 1:(obj.polinomialOrder+1)
                for j = 1:(obj.polinomialOrder+1)
                    if (i+j)<=(obj.polinomialOrder+2)
                        s = obj.computeMonomialIndeces(i,j);
                        f = matlabFunction(obj.shapeFunctions{s},'Vars',[x y]);
                        obj.shapeFunctions{s} = f;
                    end
                end
            end
        end
        
        function X = computeBasisInMonomialForm(obj)
            syms x y
            for i = 1:(obj.polinomialOrder+1)
                for j = 1:(obj.polinomialOrder+1)
                    if (i+j)<=(obj.polinomialOrder+2)
                        s = obj.computeMonomialIndeces(i,j);
                        X(s) = x^(i-1)*y^(j-1);
                    end
                end
            end
        end
        
        function B = assemblyShapeFunctionCoefficientsRHS(obj,i,j)
            I = obj.computeMonomialIndeces(i,j);
            B = zeros(obj.ndofs,1);
            B(I) = 1;
        end
        
        function A = assemblyShapeFunctionCoefficientsLHS(obj,X)
            syms x y
            A = zeros(obj.ndofs);
            for i = 1:(obj.polinomialOrder+1)
                for j = 1:(obj.polinomialOrder+1)
                    if (i+j)<=(obj.polinomialOrder+2)
                        s = obj.computeMonomialIndeces(i,j);
                        node(:) = obj.nodes(i,j,:);
                        A(s,:) = subs(X,[x y],node);
                    end
                end
            end
        end
        
        function s = computeMonomialIndeces(obj,i,j)
            n = -(obj.polinomialOrder+2);
            for m = 1:i
                n = n + obj.polinomialOrder + 1 - (m-2);
            end
            s = n+j;
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