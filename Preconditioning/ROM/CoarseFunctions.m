classdef CoarseFunctions < handle

    properties (GetAccess = public, SetAccess = private)
        f
    end

    properties (Access=private)
        mesh
        ndim
        order
        type           % "discontinous" / "continous"
    end

    methods (Access = public)

        function obj = CoarseFunctions(cParams)
            obj.init(cParams);
        end

        function fH = compute(obj)
            switch obj.ndim
                case 2
                    switch obj.type
                        case 'discontinuous'
                            fH=obj.createDiscontinuousFunction2D();
                        case 'continuous'
                            fH=obj.createContinuousFunction2D();
                    end
                case 3
                    switch obj.type
                        case 'discontinuous'
                            fH=obj.createDiscontinuous3D();
                        case 'continuous'
                            fH=obj.createContinuousFunction3D();
                    end
            end
        end

        function f=getAnalytical(obj)
            fH = obj.compute();
            nf = numel(fH);
            f  = cell(1,nf);
            switch obj.type
                case 'discontinuous'
                    nPerSegment=obj.ndim*(obj.order+1); %2fx+2fy x bmesh{i}
                    idx=1;
                    for k=1:numel(obj.mesh)
                        bMesh = obj.mesh{k}.mesh;
                        for j=1:nPerSegment
                            f{idx}=AnalyticalFunction.create(fH{idx}, bMesh);
                            idx=idx+1;
                        end
                    end

                case 'continuous'
                    bMesh=obj.mesh;
                    for i=1:nf
                        f{i}=AnalyticalFunction.create(fH{i}, bMesh);
                    end
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
                obj.ndim = obj.mesh{1}.mesh.ndim;
            else
                obj.ndim = obj.mesh.ndim;
            end

            if isfield(cParams,'type')
                obj.type = cParams.type; % continuous / discontinuous
            else
                obj.type = 'discontinuous'; % default discontinuous
            end

        end


        function f=createDiscontinuousFunction2D(obj)  
            L=obj.createBasisFunctions();
            nf = numel(obj.mesh) *obj.ndim* (obj.order + 1);
            f  = cell(1,nf);
            n  = 1;
            for k=1:numel(obj.mesh)
                bMesh=obj.mesh{k}.mesh;

                [xmax,xmin,a,b,x0,y0,~,~]=obj.NormalizeMesh(bMesh);

                if abs(xmax - xmin) < 1e-12
                    local=  @(x) (x(2,:,:)-y0)/b; %vertical
                else
                    local = @(x) (x(1,:,:)-x0)/a; %horizontal
                end
                
                for i=1:obj.order+1
                    N = @(x) L{i}( local(x) );
                    fx = @(x) [N(x); 0*x(2,:,:)];
                    fy = @(x) [0*x(1,:,:); N(x)];
                    f{1,n} = fx; n=n+1;
                    f{1,n} = fy; n=n+1;
                end
            end
        end

        function f=createContinuousFunction2D(obj)
            bMesh=obj.mesh;
            L=obj.createBasisFunctions(); 
            [~,~,a,b,x0,y0,~,~]=obj.NormalizeMesh(bMesh);
            bn=obj.getAllNodes(obj.order);

            nf=size(bn,1)*obj.ndim;
            f = cell(1,nf);
            n=1;
            for k = 1:size(bn,1)
                i = bn(k,1);
                j = bn(k,2);
                N = @(x) L{i}((x(1,:,:)-x0)/a) .* L{j}((x(2,:,:)-y0)/b);
                fx = @(x) [N(x); 0*x(2,:,:)];
                fy = @(x) [0*x(1,:,:); N(x)];
                f{1,n} = fx;   n=n+1;
                f{1,n} = fy;   n=n+1;
            end
        end


        function L=createBasisFunctions(obj)
            xi = linspace(-1, 1, obj.order+1);
            L = cell(obj.order+1, 1);
            for i = 1:obj.order+1
                L{i} = @(s) obj.lagrangeBasis(s, xi, i);
            end
        end


        function f=createDiscontinuous3D(obj)
            L=obj.createBasisFunctions();

            nf=0;
            for k=1:numel(obj.mesh)
                nf = nf + (obj.order+1)^2 * obj.ndim;
            end
            f=cell(1,nf);

            counter=1;

            for k=1:numel(obj.mesh)
                bMesh=obj.mesh{k}.mesh;
                [dir1,dir2]=obj.getFaceCoordinates(bMesh);

                for i=1:obj.order+1
                    for j=1:obj.order+1
                        N=@(x) L{i}(dir1(x)).*L{j}(dir2(x));

                        fx=@(x)[N(x);0*x(1,:,:);0*x(1,:,:)];
                        fy=@(x)[0*x(1,:,:);N(x);0*x(1,:,:)];
                        fz=@(x)[0*x(1,:,:);0*x(1,:,:);N(x)];

                        f{counter}=fx;  f{counter+1}=fy;   f{counter+2}=fz;
                        counter=counter+3;
                    end
                end

            end


        end

        function f = createContinuousFunction3D(obj)
            bMesh = obj.mesh;
            L = obj.createBasisFunctions();
            [~,~,a,b,x0,y0,c,z0] = obj.NormalizeMesh(bMesh);

            bn = obj.getAllNodes3D(obj.order);
            nf = size(bn,1) * obj.ndim;
            f = cell(1,nf);
       
            n = 1;
            for m = 1:size(bn,1)
                i = bn(m,1);
                j = bn(m,2);
                k = bn(m,3);
                N = @(x) ...
                    L{i}((x(1,:,:)-x0)/a) .* ...
                    L{j}((x(2,:,:)-y0)/b) .* ...
                    L{k}((x(3,:,:)-z0)/c);

                fx = @(x) [N(x); 0*x(1,:,:); 0*x(1,:,:)];
                fy = @(x) [0*x(1,:,:); N(x); 0*x(1,:,:)];
                fz = @(x) [0*x(1,:,:); 0*x(1,:,:); N(x)];

                f{n} = fx; n=n+1;
                f{n} = fy; n=n+1;
                f{n} = fz; n=n+1;
            end

        end

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


        function [xmax,xmin,a,b,x0,y0,c,z0]=NormalizeMesh(bMesh)
            xmax = max(bMesh.coord(:,1));
            ymax = max(bMesh.coord(:,2));
            xmin = min(bMesh.coord(:,1));
            ymin = min(bMesh.coord(:,2));
            a = (xmax - xmin)/2;
            b = (ymax - ymin)/2;
            x0 = xmin + a;
            y0 = ymin + b;

            if bMesh.ndim==3
                zmax = max(bMesh.coord(:,3));
                zmin = min(bMesh.coord(:,3));
                c = (zmax - zmin)/2;
                z0 = zmin + c;
            else
                c=[]; z0=[];
            end
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


        function bn= getAllNodes(order)
            k = order + 1;
            bn = [];

            % bottom
            for i = 1:k
                bn = [bn; i 1];
            end

            % right
            for j = 2:k
                bn = [bn; k j];
            end

            % top
            for i = k-1:-1:1
                bn = [bn; i k];
            end

            % left
            for j = k-1:-1:2
                bn = [bn; 1 j];
            end

            % interior nodes
            % for j = 2:k-1
            %     for i = 2:k-1
            %         bn = [bn; i j];
            %     end
            % end
        end

        function bn = getAllNodes3D(order)
            k = order + 1;
            bn = [];

            % 8 VERTICES
            bn = [1 1 1
                  k 1 1
                  k k 1
                  1 k 1
                  1 1 k
                  k 1 k
                  k k k
                  1 k k];

            if order == 1
                return;
            end

            t = 2:k-1;

            % EDGE NODES (12 EDGES)

            % bottom
            for i=t, bn=[bn; i 1 1]; end
            for j=t, bn=[bn; k j 1]; end
            for i=fliplr(t), bn=[bn; i k 1]; end
            for j=fliplr(t), bn=[bn; 1 j 1]; end

            % top
            for i=t, bn=[bn; i 1 k]; end
            for j=t, bn=[bn; k j k]; end
            for i=fliplr(t), bn=[bn; i k k]; end
            for j=fliplr(t), bn=[bn; 1 j k]; end

            % vertical
            for z=t, bn=[bn; 1 1 z]; end
            for z=t, bn=[bn; k 1 z]; end
            for z=t, bn=[bn; k k z]; end
            for z=t, bn=[bn; 1 k z]; end


            % FACE INTERIOR NODES
            % z=1
            for j=t
                for i=t
                    bn=[bn; i j 1];
                end
            end

            % z=k
            for j=t
                for i=t
                    bn=[bn; i j k];
                end
            end

            % y=1
            for z=t
                for i=t
                    bn=[bn; i 1 z];
                end
            end

            % y=k
            for z=t
                for i=t
                    bn=[bn; i k z];
                end
            end

            % x=1
            for z=t
                for j=t
                    bn=[bn; 1 j z];
                end
            end

            % x=k
            for z=t
                for j=t
                    bn=[bn; k j z];
                end
            end

            % VOLUME INTERIOR
            % for z=t
            %     for j=t
            %         for i=t
            %             bn=[bn; i j z];
            %         end
            %     end
            % end
        end
 
    end

end

