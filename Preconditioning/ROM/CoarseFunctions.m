classdef CoarseFunctions < handle

    properties (GetAccess = public, SetAccess = private)
        f
    end

    properties (Access=private)
        mesh
        dim
        order
        type  % "line" / "quad"
    end

    methods (Access = public)

        function obj = CoarseFunctions(cParams)
            obj.init(cParams);
        end

        function f= compute(obj)
            switch obj.type
                case 'line'
                    f=obj.createLineFunction();
                case 'quad'
                    f=obj.createQuadFunction();
                otherwise
                    error('Unknown coarse function type');
            end
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.type = cParams.type;

            if isfield(cParams,'order')
                obj.order= cParams.order;
            else
                obj.order = 1;
            end

            if iscell(obj.mesh)
                obj.dim = obj.mesh{1}.mesh.ndim;
            else
                obj.dim = obj.mesh.ndim;
            end

        end


        function f=createLineFunction(obj)  
            L=obj.createBasisFunctions();
            nf = numel(obj.mesh) *obj.dim* (obj.order + 1);
            f  = cell(1,nf);
            n  = 1;
            for k=1:numel(obj.mesh)
                bMesh=obj.mesh{k}.mesh;

                [xmax,xmin,a,b,x0,y0]=obj.NormalizeMesh(bMesh);

                if abs(xmax - xmin) < 1e-12
                    local=  @(x) (x(2,:,:)-y0)/b; %vertical
                else
                    local = @(x) (x(1,:,:)-x0)/a; %horizontal
                end
                
                for i=1:obj.order+1
                    N = @(x) L{i}( local(x) );
                    fx = @(x) [N(x); 0*x(2,:,:)];
                    fy = @(x) [0*x(1,:,:); N(x)];
                    f{1,n} = AnalyticalFunction.create(fx, bMesh); n=n+1;
                    f{1,n} = AnalyticalFunction.create(fy, bMesh); n=n+1;
                end
            end
        end

        function f=createQuadFunction(obj)
            bMesh=obj.mesh;
            L=obj.createBasisFunctions(); 
            [~,~,a,b,x0,y0]=obj.NormalizeMesh(bMesh);
            bn=obj.getBoundaryNodes(obj.order);

            nf=size(bn,1)*obj.dim;
            f = cell(1,nf);
            n=1;
            for k = 1:size(bn,1)
                i = bn(k,1);
                j = bn(k,2);
                N = @(x) L{i}((x(1,:,:)-x0)/a) .* L{j}((x(2,:,:)-y0)/b);
                fx = @(x) [N(x); 0*x(2,:,:)];
                fy = @(x) [0*x(1,:,:); N(x)];
                f{1,n} = AnalyticalFunction.create(fx, bMesh); n=n+1;
                f{1,n} = AnalyticalFunction.create(fy, bMesh); n=n+1;
            end
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


        function [xmax,xmin,a,b,x0,y0]=NormalizeMesh(bMesh)
            xmax = max(bMesh.coord(:,1));
            ymax = max(bMesh.coord(:,2));
            xmin = min(bMesh.coord(:,1));
            ymin = min(bMesh.coord(:,2));
            a = (xmax - xmin)/2;
            b = (ymax - ymin)/2;
            x0 = xmin + a;
            y0 = ymin + b;
        end

        function bn=getBoundaryNodes(order)
            k = order + 1;
            bn = zeros(4*order,2);
            n= 1;
            
            for i = 1:k                  % bottom edge
                bn(n,:) = [i, 1]; n=n+1;
            end
            for j = 2:k                  % right edge
                bn(n,:) = [k, j]; n=n+1;
            end
            for i = k-1:-1:1             % top edge
                bn(n,:) = [i, k]; n=n+1;
            end
            for j = k-1:-1:2             % left edge
                bn(n,:) = [1, j]; n=n+1;
            end

        end
    end

end

