classdef CoarseFunction < handle

    properties (GetAccess = public, SetAccess = private)
        f
    end

    methods (Access = public)

        function obj = CoarseFunction(bMesh,order)
            % order: polynomial order (1 = linear, 2 = quadratic, 3 = cubic)
            if nargin < 2
                order = 1;
            end

            % Mesh normalization
            xmax = max(bMesh.coord(:,1));
            ymax = max(bMesh.coord(:,2));
            xmin = min(bMesh.coord(:,1));
            ymin = min(bMesh.coord(:,2));
            a = (xmax - xmin)/2;
            b = (ymax - ymin)/2;
            x0 = xmin + a;
            y0 = ymin + b;
            % 1D reference node positions
            xi = linspace(-1, 1, order+1);
            L = cell(order+1, 1);
            for i = 1:order+1
                L{i} = @(s) obj.lagrangeBasis(s, xi, i);
            end
            %             % Loop only over boundary node indices
            %             f = {};
            %             for i = 1:order+1
            %                 for j = 1:order+1
            %                     if obj.isBoundaryNode(i, j, order+1)
            %                         N = @(x) L{i}((x(1,:,:)-x0)/a) .* L{j}((x(2,:,:)-y0)/b);
            %                         fx = @(x) [N(x); 0*x(2,:,:)];
            %                         fy = @(x) [0*x(1,:,:); N(x)];
            %                         f{end+1} = fx;
            %                         f{end+1} = fy;
            %                     end
            %                 end
            %             end
 % --- Detect if it's a line or full boundary ---
            if (abs(xmax - xmin) < 1e-12) || (abs(ymax - ymin) < 1e-12)
                % ✅ 1D case (single line)
                f = {};
                if abs(xmax - xmin) < 1e-12
                    % vertical line (x constant)
                    for j = 1:order+1
                        N = @(x) L{j}((x(2,:,:)-y0)/b);
                        fx = @(x) [N(x); 0*x(2,:,:)];
                        fy = @(x) [0*x(1,:,:); N(x)];
                        f{end+1} = fx;
                        f{end+1} = fy;
                    end
                else
                    % horizontal line (y constant)
                    for i = 1:order+1
                        N = @(x) L{i}((x(1,:,:)-x0)/a);
                        fx = @(x) [N(x); 0*x(2,:,:)];
                        fy = @(x) [0*x(1,:,:); N(x)];
                        f{end+1} = fx;
                        f{end+1} = fy;
                    end
                end
            else
                % ✅ 2D case (4 edges)
                n = order + 1;
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

            nfun = length(f);
            for k = 1:nfun
                uD{k} = AnalyticalFunction.create(f{k}, bMesh);
%                                 uD{k}.plot()
            end
            obj.f = uD;
        end

    end

    methods (Access = private)

        function val = lagrangeBasis(obj,s, nodes, i)
            n = length(nodes);
            val = ones(size(s));
            for j = 1:n
                if j ~= i
                    val = val .* (s - nodes(j)) / (nodes(i) - nodes(j));
                end
            end
        end

        function isB = isBoundaryNode(obj,i, j, n)
            % n = order + 1
            isB = (i == 1 || i == n || j == 1 || j == n);
        end

    end

end

