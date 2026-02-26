classdef CoarseFunction2 < handle
    properties (GetAccess = public, SetAccess = private)
        f   % cell array of AnalyticalFunction objects
    end

    methods (Access = public)
        function obj = CoarseFunction2(bMesh, order)
            if nargin < 2
                order = 1;
            end

            % --- Geometry limits ---
            xmax = max(bMesh.coord(:,1));
            ymax = max(bMesh.coord(:,2));
            xmin = min(bMesh.coord(:,1));
            ymin = min(bMesh.coord(:,2));
            a = (xmax - xmin)/2;
            b = (ymax - ymin)/2;
            x0 = xmin + a;
            y0 = ymin + b;

            % --- 1D Lagrange basis nodes on [-1,1] ---
            xi = linspace(-1,1,order+1);
            L = cell(order+1,1);
            for i = 1:order+1
                L{i} = @(s) obj.lagrangeBasis(s, xi, i);
            end

            % --- detect geometry type ---
            isVertical   = abs(xmax - xmin) < 1e-12;
            isHorizontal = abs(ymax - ymin) < 1e-12;
            isLine = isVertical || isHorizontal;

            % --- helper for zero vector field ---
            zeroF = @(x) [0*x(1,:,:); 0*x(2,:,:)];

            f = {};

            if ~isLine
                % ==========================================================
                % ✅ 2D CASE: rectangular boundary of order p
                % ==========================================================
                n = order + 1;

                % bottom edge (η = -1)
                for i = 1:n
                    N = @(x) L{i}((x(1,:,:)-x0)/a).*L{1}((x(2,:,:)-y0)/b);
                    f{end+1} = @(x) [N(x);0*x(2,:,:)];  % fx
                    f{end+1} = @(x) [0*x(1,:,:);N(x)];  % fy
                end

                % right edge (ξ = +1)
                for j = 2:n
                    N = @(x) L{n}((x(1,:,:)-x0)/a).*L{j}((x(2,:,:)-y0)/b);
                    f{end+1} = @(x) [N(x);0*x(2,:,:)]; 
                    f{end+1} = @(x) [0*x(1,:,:);N(x)];
                end

                % top edge (η = +1)
                for i = n-1:-1:1
                    N = @(x) L{i}((x(1,:,:)-x0)/a).*L{n}((x(2,:,:)-y0)/b);
                    f{end+1} = @(x) [N(x);0*x(2,:,:)]; 
                    f{end+1} = @(x) [0*x(1,:,:);N(x)];
                end

                % left edge (ξ = -1)
                for j = n-1:-1:2
                    N = @(x) L{1}((x(1,:,:)-x0)/a).*L{j}((x(2,:,:)-y0)/b);
                    f{end+1} = @(x) [N(x);0*x(2,:,:)]; 
                    f{end+1} = @(x) [0*x(1,:,:);N(x)];
                end

            else
                % ==========================================================
                % ✅ 1D CASE: horizontal or vertical line
                % ==========================================================
                n = order + 1;
                if isHorizontal
                    % vary in x, constant y
                    for i = 1:n
                        N = @(x) L{i}((x(1,:,:)-x0)/a);
                        f{end+1} = @(x) [N(x);0*x(2,:,:)]; % active edge
                        f{end+1} = zeroF;                   % fy zero
                    end
                else
                    % vary in y, constant x
                    for j = 1:n
                        N = @(x) L{j}((x(2,:,:)-y0)/b);
                        f{end+1} = zeroF;                   % fx zero
                        f{end+1} = @(x) [0*x(1,:,:);N(x)];  % active edge
                    end
                end
            end

            % --- Wrap into AnalyticalFunction objects ---
            for k = 1:length(f)
                uD{k} = AnalyticalFunction.create(f{k}, bMesh);
            end
            obj.f = uD;
        end
    end

    % ---------------------------------------------------------------------
    methods (Access = private)
        function val = lagrangeBasis(~, s, nodes, i)
            n = numel(nodes);
            val = ones(size(s));
            for j = 1:n
                if j ~= i
                    val = val .* (s - nodes(j)) / (nodes(i) - nodes(j));
                end
            end
        end
    end
end
