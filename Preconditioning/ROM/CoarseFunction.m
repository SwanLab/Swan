classdef CoarseFunction < handle

    properties (GetAccess = public, SetAccess = private)
        f % Cell array of AnalyticalFunction objects
    end

    methods (Access = public)

        function obj = CoarseFunction(bMesh, order)
            % order: polynomial order (1 = linear, 2 = quadratic, 3 = cubic)
            if nargin < 2
                order = 1;
            end

            % 1. Mesh Bounding Box
            xmax = max(bMesh.coord(:,1));
            ymax = max(bMesh.coord(:,2));
            xmin = min(bMesh.coord(:,1));
            ymin = min(bMesh.coord(:,2));

            % Calculate half-widths (handle 0 width/height later to avoid Inf)
            a = (xmax - xmin)/2;
            b = (ymax - ymin)/2;
            x0 = xmin + a;
            y0 = ymin + b;

            % 2. Precompute 1D Lagrange Basis functions
            % We precompute these to avoid re-generating them in loops
            xi = linspace(-1, 1, order+1);
            L = cell(order+1, 1);
            for k = 1:order+1
                % Create a local scope for k to ensure value capture
                L{k} = obj.createLagrange(xi, k);
            end

            raw_f = {}; % Temporary storage for function handles

            % 3. Geometry Detection
            isVertical   = abs(xmax - xmin) < 1e-12;
            isHorizontal = abs(ymax - ymin) < 1e-12;

            if isVertical || isHorizontal
                % --- 1D Case (Line) ---
                if isVertical
                    % Vertical line: x is constant, vary y (use b)
                    for j = 1:order+1
                        % Closure Fix: Pass specific L{j}, y0, b to helper
                        N = obj.create1DShape(L{j}, y0, b, 2); 
                        raw_f{end+1} = @(x) [N(x); 0*x(2,:,:)];
                        raw_f{end+1} = @(x) [0*x(1,:,:); N(x)];
                    end
                else
                    % Horizontal line: y is constant, vary x (use a)
                    for i = 1:order+1
                        N = obj.create1DShape(L{i}, x0, a, 1);
                        raw_f{end+1} = @(x) [N(x); 0*x(2,:,:)];
                        raw_f{end+1} = @(x) [0*x(1,:,:); N(x)];
                    end
                end
            else
                % --- 2D Case (Rectangle) ---
                n = order + 1;
                boundaryNodes = [];

                % Define boundary indices (CCW or similar order)
                % Bottom edge
                boundaryNodes = [boundaryNodes; [(1:n)', ones(n,1)]]; 
                % Right edge (exclude bottom node)
                boundaryNodes = [boundaryNodes; [n*ones(n-1,1), (2:n)']];
                % Top edge (exclude right node)
                boundaryNodes = [boundaryNodes; [(n-1:-1:1)', n*ones(n-1,1)]];
                % Left edge (exclude top and bottom nodes)
                boundaryNodes = [boundaryNodes; [ones(n-2,1), (n-1:-1:2)']];

                % Remove duplicates if corners were added multiple times (logic above handles it, 
                % but standard `unique` is safer if you change logic later)
                boundaryNodes = unique(boundaryNodes, 'rows', 'stable');

                for k = 1:size(boundaryNodes, 1)
                    i = boundaryNodes(k, 1);
                    j = boundaryNodes(k, 2);

                    % Closure Fix: Call helper to freeze i and j values
                    N = obj.create2DShape(L{i}, L{j}, x0, y0, a, b);

                    raw_f{end+1} = @(x) [N(x); 0*x(2,:,:)];
                    raw_f{end+1} = @(x) [0*x(1,:,:); N(x)];
                end
            end

            % 4. Wrap in AnalyticalFunction
            % Assuming AnalyticalFunction class exists
            nfun = length(raw_f);
            obj.f = cell(1, nfun);
            for k = 1:nfun
                obj.f{k} = AnalyticalFunction.create(raw_f{k}, bMesh);
            end
        end

    end

    methods (Access = private)

        % --- Factory Methods to Fix Closure Bugs ---

        function fun = createLagrange(~, nodes, i)
            % Returns the function handle for the i-th Lagrange polynomial
            fun = @(s) computeLagrange(s, nodes, i);
        end

        function fun = create1DShape(~, Li, center, scale, dim)
            % Returns N(x) for 1D. dim=1 for x-var, dim=2 for y-var.
            fun = @(x) Li( (x(dim,:,:) - center) / scale );
        end

        function fun = create2DShape(~, Li, Lj, x0, y0, a, b)
            % Returns N(x) for 2D tensor product
            fun = @(x) Li((x(1,:,:)-x0)/a) .* Lj((x(2,:,:)-y0)/b);
        end

    end
end

% --- Non-method Helper Function for cleaner Lagrange logic ---
function val = computeLagrange(s, nodes, i)
    n = length(nodes);
    val = ones(size(s));
    for j = 1:n
        if j ~= i
            val = val .* (s - nodes(j)) / (nodes(i) - nodes(j));
        end
    end
end