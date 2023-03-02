classdef CrouzeixRaviart2D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        edgeVectors
        midPoints
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = CrouzeixRaviart2D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            set(groot,'defaultAxesTickLabelInterpreter','latex');
            
            s.coord = obj.vertices;
            s.connec = [1 2 3];
            m = Mesh(s);
            
%             obj.fig = figure();
            for i=1:obj.n_vertices
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
       
        function init(obj)
            obj.computeVertices()
            obj.computeEdgesVectors()
            obj.computeMidpoints()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeEdgesVectors(obj)
            obj.edgeVectors(1,:) = obj.vertices(3,:)-obj.vertices(2,:);
            obj.edgeVectors(2,:) = obj.vertices(1,:)-obj.vertices(3,:);
            obj.edgeVectors(3,:) = obj.vertices(2,:)-obj.vertices(1,:);
        end
        
        function computeMidpoints(obj)
            for i = 1:obj.n_vertices
                if i==obj.n_vertices
                    obj.midPoints(i,:) = obj.vertices(1,:)+obj.edgeVectors(i,:)/2;
                else
                    obj.midPoints(i,:) = obj.vertices(i+1,:)+obj.edgeVectors(i,:)/2;
                end
            end
        end
        
        function computeShapeFunctions(obj)
            syms x y
            matrixLHS = [1 obj.midPoints(1,:) obj.midPoints(1,:); 1 obj.midPoints(2,:) obj.midPoints(2,:); 1 obj.midPoints(3,:) obj.midPoints(3,:)];
            for i = 1:obj.n_vertices
                vectorRHS = zeros(obj.n_vertices,1);
                vectorRHS(i) = 1;
                coefShapeFunc = vectorRHS\matrixLHS;
                obj.shapeFunctions{i} = matlabFunction(coefShapeFunc(1)+coefShapeFunc(2)*x+coefShapeFunc(3)*y,'Vars',[x y]);
            end
        end
        
    end
    
end