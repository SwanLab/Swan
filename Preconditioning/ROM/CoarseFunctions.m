classdef CoarseFunctions < handle

    properties (GetAccess = public, SetAccess = private)
        f
    end

    properties (Access=private)
        mesh
        dim
        order
        range
        activeDims
        dB
    end

    methods (Access = public)

        function obj = CoarseFunctions(cParams)
            obj.init(cParams);
        end

        function f= compute(obj)
            % order: polynomial order (1 = linear, 2 = quadratic, 3 = cubic)
            bMesh=obj.mesh;
             % Extract geometry from mesh 
            xmax = max(bMesh.coord(:,1));
            ymax = max(bMesh.coord(:,2));
            xmin = min(bMesh.coord(:,1));
            ymin = min(bMesh.coord(:,2));
            a = (xmax - xmin)/2;
            b = (ymax - ymin)/2;
            x0 = xmin + a;
            y0 = ymin + b;

            L=obj.createBasisFunctions(); 

 % --- Detect if it's a line or full boundary ---
            if (abs(xmax - xmin) < 1e-12) || (abs(ymax - ymin) < 1e-12)
                % 1D case (single line)
                f = {};
                if abs(xmax - xmin) < 1e-12
                    % vertical line (x constant)
                    for j = 1:obj.order+1
                        N = @(x) L{j}((x(2,:,:)-y0)/b);
                        fx = @(x) [N(x); 0*x(2,:,:)];
                        fy = @(x) [0*x(1,:,:); N(x)];
                        f{end+1} = fx;
                        f{end+1} = fy;
                    end
                else
                    % horizontal line (y constant)
                    for i = 1:obj.order+1
                        N = @(x) L{i}((x(1,:,:)-x0)/a);
                        fx = @(x) [N(x); 0*x(2,:,:)];
                        fy = @(x) [0*x(1,:,:); N(x)];
                        f{end+1} = fx;
                        f{end+1} = fy;
                    end
                end
            else
                % 2D case (4 edges)
                n = obj.order + 1;
                boundaryNodes = [];

                % bottom edge
                for i = 1:n
                    boundaryNodes(end+1,:) = [i, 1];
                end
                % right edge
                for j = 2:n
                    boundaryNodes(end+1,:) = [n, j];
                end
                % top edge
                for i = n-1:-1:1
                    boundaryNodes(end+1,:) = [i, n];
                end
                % left edge
                for j = n-1:-1:2
                    boundaryNodes(end+1,:) = [1, j];
                end

                f = {};
                for k = 1:size(boundaryNodes,1)
                    i = boundaryNodes(k,1);
                    j = boundaryNodes(k,2);
                    N = @(x) L{i}((x(1,:,:)-x0)/a) .* L{j}((x(2,:,:)-y0)/b);
                    fx = @(x) [N(x); 0*x(2,:,:)];
                    fy = @(x) [0*x(1,:,:); N(x)];
                    f{end+1} = fx;
                    f{end+1} = fy;
                end
            end

            uD= cell(size(f));
            for k = 1:length(f)
                uD{k} = AnalyticalFunction.create(f{k}, bMesh);
%                                 uD{k}.plot()
            end
            f = uD;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh =cParams.mesh;
            if isfield(cParams,'order')
                obj.order= cParams.order;
            else
                obj.order = 1;
            end

            obj.dim= size(cParams.mesh.coord,2);
            obj.range = max(obj.mesh.coord) - min(obj.mesh.coord);
            obj.activeDims = find(obj.range > 1e-12);
            obj.dB = numel(obj.activeDims);
        end

        function L=createBasisFunctions(obj)
            xi = linspace(-1, 1, obj.order+1);
            L = cell(obj.order+1, 1);
            for i = 1:obj.order+1
                L{i} = @(s) obj.lagrangeBasis(s, xi, i);
            end
        end


    end

    methods (Access=private,Static)
        function val = lagrangeBasis(s, nodes, i)
            n = length(nodes);
            val = ones(size(s));
            for j = 1:n
                if j ~= i
                    val = val .* (s - nodes(j)) / (nodes(i) - nodes(j));
                end
            end
        end
    end

end

