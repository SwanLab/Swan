classdef LagrangeSimplicial2D < handle
   
    properties (Access = private)
        k
    end
    
    
    properties (Access = public)
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
            obj.fig = figure();
            syms x y
            hold on
            if obj.k==1
                for i=1:obj.n_nodes
                    subplot(1,3,i)
                    func = piecewise(x+y<=1,obj.shapeFunctions{i});
                    fsurf(func,[0 1])
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    title("i:"+obj.nodes.index(i))
                    grid on
                end
            elseif obj.k == 2
                for i=1:obj.n_nodes
                    subplot(2,3,i)
                    func = piecewise(x+y<=1,obj.shapeFunctions{i});
                    fsurf(func,[0 1])
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    title("i:"+obj.nodes.index(i,1)+", j:"+obj.nodes.index(i,2))
                    grid on
                end
            elseif obj.k == 3
                for i=1:obj.n_nodes
                    subplot(2,5,i)
                    func = piecewise(x+y<=1,obj.shapeFunctions{i});
                    fsurf(func,[0 1])
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    zlim([-0.5 1])
                    title("i:"+obj.nodes.index(i,1)+", j:"+obj.nodes.index(i,2)+", k:"+obj.nodes.index(i,3))
                    grid on
                end
            end
            hold off
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj,k)
            obj.k = k;
            obj.computeVertices()
            obj.computeNormalVectors()
            obj.computeNodes()
            obj.computeBarycentricCoords()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = [sqrt(2)/2,sqrt(2)/2;0,-1;-1,0];
        end
        
        function computeNodes(obj)
            obj.n_nodes = nchoosek(2+obj.k,obj.k);
            
            if obj.k == 1
                obj.nodes.coord = [0,0;0,1;1,0];
                obj.nodes.index = [1,2,3];
            elseif obj.k == 2
                obj.nodes.coord = [0,0;0,1;0,0.5;0.5,0;0.5,0.5;1,0];
                obj.nodes.index = [1,0;1,2;2,0;1,3;2,3;3,0];
            elseif obj.k == 3
                obj.nodes.coord = [0,0;0,1/3;0,2/3;0,1;1/3,0;1/3,1/3;1/3,2/3;2/3,0;2/3,1/3;2/3,1/3;1,0];
                obj.nodes.index = [1,0,0;1,2,0;2,1,0;2,0,0;1,3,0;1,2,3;2,3,0;3,1,0;3,2,0;3,0,0];
            end
        end
        
        function computeBarycentricCoords(obj)
            syms x y
            for i = 1:obj.n_vertices
                if i~=obj.n_vertices
                    obj.barycentricCoords{i} = symfun(1-dot([x,y]-obj.vertices(i,:),obj.normalVectors(i,:))/dot(obj.vertices(i+1,:)-obj.vertices(i,:),obj.normalVectors(i,:)),[x,y]);
                else
                    obj.barycentricCoords{i} = symfun(1-dot([x,y]-obj.vertices(i,:),obj.normalVectors(i,:))/dot(obj.vertices(1,:)-obj.vertices(i,:),obj.normalVectors(i,:)),[x,y]);
                end
            end
        end
        
        function computeShapeFunctions(obj)
            if obj.k==1
                obj.shapeFunctions = obj.barycentricCoords;
                
            elseif obj.k==2
                obj.shapeFunctions = cell(obj.n_nodes,1);
                for i = 1:obj.n_nodes
                    if obj.nodes.index(i,2) == 0
                        obj.shapeFunctions{i} = obj.barycentricCoords{obj.nodes.index(i,1)}*(2*obj.barycentricCoords{obj.nodes.index(i,1)}-1);
                    else
                        obj.shapeFunctions{i} = 4*obj.barycentricCoords{obj.nodes.index(i,1)}*obj.barycentricCoords{obj.nodes.index(i,2)};
                    end
                end
                
            elseif obj.k==3
                obj.shapeFunctions = cell(obj.n_nodes,1);
                for i = 1:obj.n_nodes
                    if obj.nodes.index(i,2) == 0 && obj.nodes.index(i,3) == 0
                        func = obj.barycentricCoords{obj.nodes.index(i,1)};
                        obj.shapeFunctions{i} = 0.5*func*(3*func-1)*(3*func-2);
                    elseif obj.nodes.index(i,2) ~= 0 && obj.nodes.index(i,3) == 0
                        func_i = obj.barycentricCoords{obj.nodes.index(i,1)};
                        func_j = obj.barycentricCoords{obj.nodes.index(i,2)};
                        obj.shapeFunctions{i} = 9/2*func_i*(3*func_i-1)*func_j;
                    else
                        func_i = obj.barycentricCoords{obj.nodes.index(i,1)};
                        func_j = obj.barycentricCoords{obj.nodes.index(i,2)};
                        func_k = obj.barycentricCoords{obj.nodes.index(i,3)};
                        obj.shapeFunctions{i} = 27*func_i*func_j*func_k;
                    end
                end
            end
        end
        
    end
    
end