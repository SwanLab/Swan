classdef CrouzeixRaviart1D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        ndofs
        shapeFunctions
        fig
        simplicial
    end
    
    
    methods (Access = public)
    
        function obj = CrouzeixRaviart1D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
%             obj.fig = figure();
            for i=1:obj.n_vertices
%                 subplot(1,2,i)
                figure();
                fplot(obj.shapeFunctions{i},[0 1]);
                xlabel('x'); ylabel('y');
                title("i: "+string(i-1));
                grid on
            end
            hold off
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.simplicial = Simplicial1D();
            obj.computeVertices();
            obj.computeNdof();
            obj.computeShapeFunctions();
        end
        
        function computeVertices(obj)
            obj.vertices = obj.simplicial.vertices;
            obj.n_vertices = length(obj.vertices);
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.vertices);
        end
        
        function computeShapeFunctions(obj)
            syms x
            ndof = obj.ndofs;
            shapeFunc = cell(ndof,1);
            for i = 1:ndof
                LHS = obj.applyLinearForm();
                RHS = obj.computeLinearFormValues(i);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                shapeFunc{i} = matlabFunction(c(1)+c(2)*x);
            end
            obj.shapeFunctions = shapeFunc;
        end
        
        function c = computeShapeFunctionCoefficients(~,LHS,RHS)
            c = LHS\RHS;
        end
        
        function LHS = applyLinearForm(obj)
            LHS = [1 obj.vertices(2); 1 obj.vertices(1)];
        end
        
        function RHS = computeLinearFormValues(obj,i)
            RHS = zeros(obj.n_vertices,1);
            RHS(i) = 1;
        end
        
    end
    
end