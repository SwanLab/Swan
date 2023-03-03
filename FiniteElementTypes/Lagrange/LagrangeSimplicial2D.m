classdef LagrangeSimplicial2D < handle
   
    properties (Access = public)
        polynomialOrder
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
            obj.polynomialOrder = k;
            obj.computeVertices();
            obj.computeNdof();
            obj.computeNodes();
            obj.computeShapeFunctions();
        end

        function computeNdof(obj)
            k = obj.polynomialOrder;                                    
            n = nchoosek(2+k,k);            
            obj.ndofs = n;
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeNodes(obj)
            k = obj.polynomialOrder;                                    
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
            k = obj.polynomialOrder;
            X = obj.computeBasisInMonomialForm();
            shapeFuncs = cell(obj.ndofs,1);
            for i = 1:(k+1)
                for j = 1:(k+1)
                    if (i+j)<=(k+2)
                        a = obj.computeShapeFunctionCoefficients(X,i,j);
                        s = obj.computeMonomialIndeces(i,j);
                        shapeFuncs{s} = matlabFunction(X*a,'Vars',[x y]);
                    end
                end
            end
            
            obj.shapeFunctions = shapeFuncs;
        end
        
        function s = computeShapeFunctionCoefficients(obj,X,i,j)
            A = obj.applyLinearFormInMonomialForm(X);
            b = obj.computeLinearFormValues(i,j);
            s = A\b;
        end
        
        function X = computeBasisInMonomialForm(obj)
            k = obj.polynomialOrder;
            syms x y
            for i = 1:(k+1)
                for j = 1:(k+1)
                    if (i+j)<=(k+2)
                        s = obj.computeMonomialIndeces(i,j);
                        X(s) = x^(i-1)*y^(j-1);
                    end
                end
            end
        end
        
        function B = computeLinearFormValues(obj,i,j)
            I = obj.computeMonomialIndeces(i,j);
            B = zeros(obj.ndofs,1);
            B(I) = 1;
        end
        
        function A = applyLinearFormInMonomialForm(obj,X)
            k = obj.polynomialOrder;
            syms x y
            A = zeros(obj.ndofs);
            for i = 1:(k+1)
                for j = 1:(k+1)
                    if (i+j)<=(k+2)
                        s = obj.computeMonomialIndeces(i,j);
                        node(:) = obj.nodes(i,j,:);
                        A(s,:) = subs(X,[x y],node);
                    end
                end
            end
        end
        
        function s = computeMonomialIndeces(obj,i,j)
            k = obj.polynomialOrder;
            n = -(k+2);
            for m = 1:i
                n = n + k + 1 - (m-2);
            end
            
            s = n+j;
        end
        
    end
    
end