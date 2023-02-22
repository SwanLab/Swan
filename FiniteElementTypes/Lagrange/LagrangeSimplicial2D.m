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
            m.plot()

            a = trisurf(m.connec,s.coord(:,1),s.coord(:,2),obj.shapeFunctions{1,1});

            obj.fig = figure();
            syms x y
            hold on
            for i=1:obj.k+1
                for j = 1:obj.k+1
                    if (j+i)<=(obj.k+2)
                        subplot(obj.k+1,obj.k+1,(i-1)*(obj.k+1)+j)
                        func = piecewise(x+y<=1,obj.shapeFunctions{i,j});
                        fsurf(func,[0 1])
                        zlim([-0.5 1])
                        xlabel('x')
                        ylabel('y')
                        zlabel('z')
                        title("i:"+string(i-1)+", j:"+string(j-1))
                        grid on
                    end
                end
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj,k)
            obj.k = k;
            obj.computeVertices()
            obj.computeNodes()
            obj.computeShapeFunctions()
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
            
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    if (i+j)<=(obj.k+2)
                        n = -(obj.k+2);
                        for m = 1:i
                            n = n + obj.k + 1 - (m-2);
                        end
                        I = n+j;
                        X(I) = x^(i-1)*y^(j-1);
                    end
                end
            end
                      
            for i = 1:(obj.k+1)
                for j = 1:(obj.k+1)
                    if (i+j)<=(obj.k+2)
                        n = -(obj.k+2);
                        for m = 1:i
                            n = n + obj.k + 1 - (m-2);
                        end
                        I = n+j;

                        b = zeros(obj.n_nodes,1);
                        b(I) = 1;
                        
                        A = zeros(obj.n_nodes);
                        for ii = 1:(obj.k+1)
                            for jj = 1:(obj.k+1)
                                if (ii+jj)<=(obj.k+2)
                                    nn = -(obj.k+2);
                                    for m = 1:ii
                                        nn = nn + obj.k + 1 - (m-2);
                                    end
                                    II = nn+jj;
                                    node(:) = obj.nodes(ii,jj,:);
                                    A(II,:) = subs(X,[x y],node);
                                end
                            end
                        end
                        
                        s = A\b;
                        obj.shapeFunctions{i,j} = X*s;
                    end
                end
            end
        end
        
    end
    
end