classdef LagrangeSimplicial1D < handle
   
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
    
        function obj = LagrangeSimplicial1D(polynomialOrder)
            obj.init(polynomialOrder);
        end
    
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            set(groot,'defaultAxesTickLabelInterpreter','latex');
            
            k = obj.polynomialOrder;
%             obj.fig = figure();
            hold on
            for i = 1:k+1
%                 subplot(1,k+1,i)
                figure();
                fplot(obj.shapeFunctions{i},[0 1]);
                xlabel('x'); ylabel('y');
                title("$i="+string(i-1)+"$");
                grid on
            end
            hold off
        end
        
    end
    
    
    methods (Access = private)
       
        function init(obj,polynomialOrder)
            obj.polynomialOrder = polynomialOrder;
            obj.computeVertices();
            obj.computeNdof();
            obj.computeNodes();
            obj.computeShapeFunctions();
        end
        
        function computeVertices(obj)
            obj.n_vertices = 2;
            obj.vertices = [0,1];
        end
        
        function computeNdof(obj)
            obj.ndofs = nchoosek(1+k,k);
        end
        
        function computeNodes(obj)
            k = obj.polynomialOrder;
            ndof = obj.ndofs;
            coord = zeros(ndof,1); index = zeros(ndof,1);
            for i=1:ndof
                coord(i) = (i-1)/k;
                index(i) = i;
            end
            
            obj.nodes.coord = coord;
            obj.nodes.index = index;
        end

        function computeShapeFunctions(obj)
            for i = 1:obj.ndofs
                ndof = obj.ndofs;
                syms x
                shapeFunc = 1;
                for j = 1:ndof
                    if i~=j
                        shapeFunc = shapeFunc * (x-obj.nodes.coord(j))/(obj.nodes.coord(i)-obj.nodes.coord(j));
                    end
                end
                
                obj.shapeFunctions{i} = matlabFunction(shapeFunc);
            end
        end
        
    end
    
end