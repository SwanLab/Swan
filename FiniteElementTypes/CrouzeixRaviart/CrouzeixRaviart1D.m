classdef CrouzeixRaviart1D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        n_nodes
        nodes
        edgeVectors
        midPoints
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = CrouzeixRaviart1D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            hold on
            for i=1:obj.n_vertices
                subplot(1,2,i)
                fplot(obj.shapeFunctions{i},[0 1])
                xlabel('x')
                ylabel('y')
                title("i:"+i)
                grid on
            end
            hold off
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.computeVertices()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 2;
            obj.vertices = [0,1];
        end
        
        function computeShapeFunctions(obj)
            syms x
            for i = 1:obj.n_vertices
                matrixLHS = asssemblyLHS();
                vectorRHS = assemblyRHS(i);
                coefShapeFunc = matrixLHS\vectorRHS;
                obj.shapeFunctions{i} = matlabFunction(coefShapeFunc(1)+coefShapeFunc(2)*x);
            end
        end
        
        function matrixLHS = assemblyLHS(obj)
            matrixLHS = [1 obj.vertices(2); 1 obj.vertices(1) ];
        end
        
        function vectorRHS = assemblyRHS(obj,i)
            vectorRHS = zeros(obj.n_vertices,1);
            vectorRHS(i) = 1;
        end
        
    end
    
end