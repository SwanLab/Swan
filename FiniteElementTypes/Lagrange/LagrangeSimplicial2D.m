classdef LagrangeSimplicial2D < handle
   
    properties (Access = private)
        xSym
        ySym
        shapeFunctionsSym
        shapeFunctionsDiffSym
        domainK
        nodes
        polynomialOrder
    end
    
    properties (Access = public)
        ndofs
        shapeFunctions
        shapeFunctionsDiff
    end
    
    
    methods (Access = public)
    
        function obj = LagrangeSimplicial2D(k)
            obj.init(k);
        end
        
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            ndof = obj.ndofs;
            
            m = obj.createPlotMesh();
            for s=1:ndof
                figure();
                m.plot();
                trisurf(m.connec,m.coord(:,1),m.coord(:,2),obj.shapeFunctions{s}(m.coord(:,1),m.coord(:,2)));
                shading FLAT
                xlim([0 1]); ylim([0 1]);
                xlabel('x'); ylabel('y'); zlabel('z');
                title("Shape function (s = "+string(s-1)+")");
                grid on
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj,k)
            obj.polynomialOrder = k;
            obj.domainK = Simplicial2D();
            
            obj.xSym = sym('x','real');
            obj.ySym = sym('y','real');
            
            obj.computeNdof();
            obj.computeNodes();
            obj.computeShapeFunctionsSym();
            obj.computeShapeFunctions();
            obj.computeShapeFunctionsDiffSym();
            obj.computeShapeFunctionsDiff();
        end

        function computeNdof(obj)
            k = obj.polynomialOrder;                                    
            n = nchoosek(2+k,k);            
            obj.ndofs = n;
        end
        
        function computeNodes(obj)
            k = obj.polynomialOrder;
            ndof = obj.ndofs;
            nod = zeros(ndof,2);
            
            for j = 1:k+1
                for i = 1:k+1
                    if (i+j)<=(k+2)
                        s = obj.computeMonomialIndeces(i,j);
                        nod(s,1)=(i-1)/k;
                        nod(s,2)=(j-1)/k;
                    end
                end
            end
            
            obj.nodes = nod;
        end
        
        function computeShapeFunctionsSym(obj)
            shFunc = cell(obj.ndofs,1);
            ndof = obj.ndofs;
            
            basisMSym = obj.computeBasisInMonomialForm();
            basisM = matlabFunction(basisMSym,'Vars',[obj.xSym obj.ySym]);
            for s = 1:ndof
                c = obj.computeShapeFunctionCoefficients(basisM,s);
                shFunc{s} = basisMSym*c;
            end
            
            obj.shapeFunctionsSym = shFunc;
        end
        
        function computeShapeFunctions(obj)
            ndof = obj.ndofs;
            shFunc = cell(ndof,1);
            shFuncSym = obj.shapeFunctionsSym;
            x = obj.xSym;
            y = obj.ySym;
            
            for s = 1:ndof
                shFunc{s} = matlabFunction(shFuncSym{s},'Vars',[x y]);
            end
            
            obj.shapeFunctions = shFunc;
        end
        
        function c = computeShapeFunctionCoefficients(obj,X,s)
            LHS = obj.applyLinearFormInMonomialForm(X);
            RHS = obj.computeLinearFormValues(s);
            c = LHS\RHS;
        end
        
        function basisMSym = computeBasisInMonomialForm(obj)
            k = obj.polynomialOrder;
            x = obj.xSym;
            y = obj.ySym;
            
            for j = 1:(k+1)
                for i = 1:(k+1)
                    if (i+j)<=(k+2)
                        s = obj.computeMonomialIndeces(i,j);
                        basisMSym(s) = x^(i-1)*y^(j-1);
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
        
        function s = computeMonomialIndeces(obj,j,i)
            k = obj.polynomialOrder;
            
            n = -(k+2);
            for m = 1:i
                n = n + k + 1 - (m-2);
            end
            
            s = n+j;
        end
        
        function m = createPlotMesh(obj)
            s.coord = obj.domainK.vertices;
            s.connec = [1 2 3];
            m = Mesh(s);
            for i=1:4
                m = m.remesh(2);
            end
        end
        
        function computeShapeFunctionsDiffSym(obj)
            f = obj.shapeFunctionsSym;
            shD = cell(length(f),2);
            ndof = obj.ndofs;
            x = obj.xSym;
            y = obj.ySym;
            
            for s = 1:ndof
                shD{s,1} = diff(f{s},x);
                shD{s,2} = diff(f{s},y);
            end
            
            obj.shapeFunctionsDiffSym = shD;
        end
        
        function computeShapeFunctionsDiff(obj)
            ndof = obj.ndofs;
            shD = cell(ndof,2);
            shDSym = obj.shapeFunctionsDiffSym;
            x = obj.xSym;
            y = obj.ySym;
            
            for s = 1:ndof
                shD{s,1} = matlabFunction(shDSym{s,1},'Vars',[x y]);
                shD{s,2} = matlabFunction(shDSym{s,2},'Vars',[x y]);
            end
            
            obj.shapeFunctionsDiff = shD;
        end
        
    end
    
end