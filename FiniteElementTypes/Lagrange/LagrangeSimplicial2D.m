classdef LagrangeSimplicial2D < handle
   
    properties (Access = public)
        polynomialOrder
        vertices
        ndofs
        nodes
        shapeFunctions
        fig
        simplicial
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
            obj.simplicial = Simplicial2D();
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
            obj.vertices = obj.simplicial.vertices;
        end
        
        function computeNodes(obj)
            k = obj.polynomialOrder;
            ndof = obj.ndofs;
            node = zeros(ndof,2);
            for i = 1:k+1
                for j = 1:k+1
                    if (i+j)<=(k+2)
                        s = obj.computeMonomialIndeces(i,j);
                        node(s,1)=(i-1)/k;
                        node(s,2)=(j-1)/k;
                    end
                end
            end
            obj.nodes = node;
        end
        
        function computeShapeFunctions(obj)
            syms x y
            basisMonomialFormSym = obj.computeBasisInMonomialForm();
            basisMonomialForm = matlabFunction(basisMonomialFormSym);
            shapeFuncs = cell(obj.ndofs,1);
            for s = 1:obj.ndofs
                c = obj.computeShapeFunctionCoefficients(basisMonomialForm,s);
                shapeFuncs{s} = matlabFunction(basisMonomialFormSym*c,'Vars',[x y]);
            end
            
            obj.shapeFunctions = shapeFuncs;
        end
        
        function coefs = computeShapeFunctionCoefficients(obj,X,s)
            LHS = obj.applyLinearFormInMonomialForm(X);
            RHS = obj.computeLinearFormValues(s);
            coefs = LHS\RHS;
        end
        
        function basisMonomialFormSym = computeBasisInMonomialForm(obj)
            k = obj.polynomialOrder;
            syms x y
            for i = 1:(k+1)
                for j = 1:(k+1)
                    if (i+j)<=(k+2)
                        s = obj.computeMonomialIndeces(i,j);
                        basisMonomialFormSym(s) = x^(i-1)*y^(j-1);
                    end
                end
            end
        end
        
        function RHS = computeLinearFormValues(obj,s)
            RHS = zeros(obj.ndofs,1);
            RHS(s) = 1;
        end
        
        function LHS = applyLinearFormInMonomialForm(obj,X)
            LHS = zeros(obj.ndofs);
            for s = 1:obj.ndofs
                node(:) = obj.nodes(s,:);
                LHS(s,:) = X(node(1),node(2));
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