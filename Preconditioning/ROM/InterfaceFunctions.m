function f = InterfaceFunctions(mesh, order)
   % Ensure default order
    if nargin < 2
        order = 1;
    end

    % --- Extract geometry from mesh ---
    xmax = max(mesh.coord(:,1));
    ymax = max(mesh.coord(:,2));
    xmin = min(mesh.coord(:,1));
    ymin = min(mesh.coord(:,2));

    a  = (xmax - xmin)/2;
    b  = (ymax - ymin)/2;
    x0 = xmin + a;
    y0 = ymin + b;

    % --- Reference 1D nodes and basis functions ---
    xi = linspace(-1, 1, order + 1);
    L  = cell(order + 1, 1);
    for i = 1:order+1
        ii = i;  % capture index to fix MATLAB closure issue
        L{ii} = @(s) lagrange1D(s, xi, ii);
    end

    % --- Construct tensor-product shape functions ---
    N = cell(order + 1, order + 1);
    for i = 1:order+1
        for j = 1:order+1
            ii = i; jj = j;
            N{ii,jj} = @(x) ...
                L{ii}((x(1,:,:)-x0)/a) .* ...
                L{jj}((x(2,:,:)-y0)/b);
        end
    end

    % --- Collect boundary nodes (counterclockwise) ---
    n = order + 1;
    boundaryNodes = [];
    % bottom (η = -1)
    for i = 1:n
        boundaryNodes(end+1,:) = [i, 1];
    end
    % right (ξ = +1)
    for j = 2:n
        boundaryNodes(end+1,:) = [n, j];
    end
    % top (η = +1)
    for i = n-1:-1:1
        boundaryNodes(end+1,:) = [i, n];
    end
    % left (ξ = -1)
    for j = n-1:-1:2
        boundaryNodes(end+1,:) = [1, j];
    end

    % --- Vector-valued FE functions ---
    f = cell(1, 2*size(boundaryNodes,1));
    for k = 1:size(boundaryNodes,1)
        i = boundaryNodes(k,1);
        j = boundaryNodes(k,2);
        Ni = N{i,j};

        fx = @(x) [Ni(x); 0*x(2,:,:)];  % x-component shape
        fy = @(x) [0*x(1,:,:); Ni(x)];  % y-component shape

        f{2*k-1} = fx;
        f{2*k}   = fy;
    end

       nfun = length(f);
   for k = 1:nfun
       uD{k} = AnalyticalFunction.create(f{k}, mesh);
   end
   f = uD;
end


function L = lagrange1D(s, nodes, i)
    % Returns i-th Lagrange polynomial evaluated at s
    n = numel(nodes);
    L = ones(size(s));
    for j = 1:n
        if j ~= i
            L = L .* (s - nodes(j)) / (nodes(i) - nodes(j));
        end
    end
end