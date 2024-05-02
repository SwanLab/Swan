classdef UnsteadyConvectionProblem < handle

    properties (Access = private)
        method
        revolutions
        timeSteps
    end

    properties (Access = private)
        nElemX
        nElemY
        coord
        connec
        a
    end

    methods (Access = public)
        function obj = UnsteadyConvectionProblem(cParams)
            obj.init(cParams);
            obj.setHardCodedProperties();
        end

        function u = compute(obj)
            nx    = obj.nElemX;
            ny    = obj.nElemY;
            X     = obj.coord;
            T     = obj.connec;
            numnp = size(X,1);
            Conv  = obj.a;
            meth  = obj.method;
            n_end = obj.revolutions;
            nstep = obj.timeSteps;
            % Number of Gauss points (numerical quadrature)
            ngaus=4;
            % Quadrature
            [pospg,wpg] = Quadrature(ngaus);
            % Shape Functions
            [N,Nxi,Neta] = ShapeFunc(pospg);

            % COMPUTATION OF THE MATRICES
            % Matrices obtained by discretizing the Galerkin weak-form
            M = CreMassMat(X,T,pospg,wpg,N,Nxi,Neta);
            C = CreConvMat(X,T,Conv,pospg,wpg,N,Nxi,Neta);
            K = CreStiffMat(X,T,Conv,pospg,wpg,N,Nxi,Neta);
            % Vectors related to source term
            v1 = Crevect1(X,T,pospg,wpg,N,Nxi,Neta);          % (w,s)
            v2 = Crevect2(X,T,Conv,pospg,wpg,N,Nxi,Neta);     % (a·grad(w),s)
            % Matrices obtained by discretizing the terms on the outflow boundary
            % Elements on the outflow boundary
            midx = round(nx/2);
            midy = round(ny/2);
            elemy0_out = 1 : midx;
            elemx1_out = nx : nx : nx*midy;
            elemy1_out = nx*ny : -1 : nx*(ny-1)+1+midx;
            elemx0_out = nx*(ny-1)+1 : -nx : nx*midy+1;
            %
            Mo = CreOutMat1 (X,T,Conv,elemy0_out,[1,2]);
            Mo = Mo + CreOutMat1 (X,T,Conv,elemx1_out,[2,3]);
            Mo = Mo + CreOutMat1 (X,T,Conv,elemy1_out,[3,4]);
            Mo = Mo + CreOutMat1 (X,T,Conv,elemx0_out,[4,1]);
            %
            Co = CreOutMat2 (X,T,Conv,elemy0_out,[1,2]);
            Co = Co + CreOutMat2 (X,T,Conv,elemx1_out,[2,3]);
            Co = Co + CreOutMat2 (X,T,Conv,elemy1_out,[3,4]);
            Co = Co + CreOutMat2 (X,T,Conv,elemx0_out,[4,1]);
            %
            vo = CreOutVect (X,T,Conv,elemy0_out,[1,2]);
            vo = vo + CreOutVect (X,T,Conv,elemx1_out,[2,3]);
            vo = vo + CreOutVect (X,T,Conv,elemy1_out,[3,4]);
            vo = vo + CreOutVect (X,T,Conv,elemx0_out,[4,1]);

            % DATA FOR THE TRANSIENT ANALYSIS
            t_end = 2*pi*n_end;
            dt = t_end/nstep;

            % INITIAL CONDITION: COSINE HILL
            u = zeros(numnp,nstep+1);
            sigma = 0.2; xref=1/6;
            for i=1:numnp
                xdim = (X(i,1)-xref)/sigma;
                ydim = (X(i,2)-xref)/sigma;
                test = xdim^2 + ydim^2;
                if test>1
                    u(i,1) = 0.0;
                else
                    u(i,1) = 0.25*(1+cos(pi*xdim))*(1+cos(pi*ydim));
                end
            end

            % COMPUTATION OF MATRICES FOR THE TRANSIENT ANALYSIS
            if meth == 1
                A = M;
                B = dt*(C - (dt/2)*K - Mo + (dt/2)*Co);
                f = dt*(v1 + (dt/2)*(v2-vo));
            elseif meth == 2
                Md = diag(M*ones(numnp,1));
                Mod = diag(Mo*ones(numnp,1));
                A = Md;
                B = dt*(C - (dt/2)*K - Mod + (dt/2)*Co);
                f = dt*(v1 + (dt/2)*(v2-vo));
            elseif meth == 3
                A = M +(dt^2/6)*(K - Co);
                B = dt*(C - (dt/2)*K - Mo + (dt/2)*Co);
                f = dt*((dt/2)*(v2 - vo) +v1);
            elseif meth == 4
                A = M - (dt/2)*C + (dt/2)*Mo;
                B = dt*C - dt*Mo;
                f = dt*v1;
            elseif meth == 5
                Md = diag(M*ones(numnp,1));
                Mod = diag(Mo*ones(numnp,1));
                A = Md - (dt/2)*C + (dt/2)*Mod;
                B = dt*C - dt*Mod;
                f = dt*v1;
            elseif meth == 6
                A = M + (dt/2)*(C + C') + (dt^2/4)*K ;
                B = -dt*(C' + (dt/2)*K);
                f = dt*v1 + (dt^2/2)*v2;
            elseif meth == 7
                alpha = 1/9;
                A1 = M;
                B1 = -(dt/3)*C'- alpha*dt^2*(K - Co);
                f1 = (dt/3)*v1 + alpha*dt^2*(v2 - vo);
                A2 = M;
                B2 = -dt*C';
                C2 = - (dt^2/2)*(K-Co);
                f2 = dt*v1 - (dt^2/2)*(v2 - vo);
            else
                error('Unavailable method')
            end


            % ESSENTIAL BOUNDARY CONDITIONS
            % Homogeneous Dirichlet boundary conditions are imposed on the inflow boundary
            % using Lagrange multipliers
            %
            % Nodes on the inflow boundary
            midx = round((nx+1)/2);
            midy = round((ny+1)/2);
            nodesDy0_in = (midx : nx+1)';
            nodesDx1_in = (midy*(nx+1) : (nx+1) : (ny+1)*(nx+1))';
            nodesDy1_in = ((ny+1)*(nx+1)- midx+1 : -1 : ny*(nx+1)+1)';
            nodesDx0_in = ((midy-1)*(nx+1)+1 : -(nx+1) : 1)';
            nodesD_in   = [nodesDy0_in; nodesDx1_in; nodesDy1_in; nodesDx0_in];
            %
            nDir = length(nodesD_in);
            CDir = [nodesD_in,zeros(nDir,1)];
            %
            Accd  = zeros(nDir,numnp);
            Accd(:,CDir(:,1)) = eye(nDir);
            bccd = CDir(:,2);
            % Total matrix (including boundary conditions)
            if meth == 7    % 2-step method
                A1tot = [A1 Accd';Accd zeros(nDir,nDir)];
                [L1,U1] = lu(A1tot);
                A2tot = [A2 Accd';Accd zeros(nDir,nDir)];
                [L2,U2] = lu(A2tot);
            else
                Atot = [A Accd';Accd zeros(nDir,nDir)];
                [L,U] = lu(Atot);
            end

            % TRANSIENT SOLUTION
            for n = 1 : nstep
                if meth == 7    % 2-step method
                    btot = [B1*u(:,n)+ f1; bccd];
                    aux =  U1\(L1\btot);
                    u_m =  u(:,n) + aux(1:numnp);
                    btot = [B2*u(:,n) + C2*u_m + f2; bccd];
                    aux  = U2\(L2\btot);
                    u(:,n+1) = u(:,n) + aux(1:numnp);
                else
                    btot = [B*u(:,n)+f; bccd];
                    aux  = U\(L\btot);
                    u(:,n+1) = u(:,n) + aux(1:numnp);
                end
            end
            if not(max(u(:,nstep+1))<100 &&  min(u(:,nstep+1))>-100)
                error('Time step too large');
            end
        end

        function plot(obj,u,i)
            x_lo=-1/2; x_up=1/2; y_lo=-1/2; y_up=1/2;
            X     = obj.coord;
            nx    = obj.nElemX;
            ny    = obj.nElemY;
            % Solution at time t=t_end
            figure(1); clf;
            set(gca,'FontSize',12,...
                'XTick', [-0.5,0,0.5],'YTick', [-0.5,0,0.5],'ZTick', 0:0.25:1);
            [xx,yy,sol] = MatSol(X,nx,ny,u(:,i));
            surface(xx,yy,sol);
            view([40,30])
            axis([x_lo,x_up,y_lo,y_up,0,1])
            grid on;

            % Contour plot of the solution at time t = t_end
            figure(2); clf;
            set(gca,'FontSize',12);
            [C,h]=contour(xx,yy,sol,[-0.1,-0.01,0.1:0.1:1.0]);
            clabel(C,h);
            axis equal; axis([x_lo,x_up,y_lo,y_up]);
            set(gca,'XTick',[x_lo,(x_up+x_lo)/2,x_up]);
            set(gca,'YTick',[y_lo,(y_up+y_lo)/2,y_up]);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.method      = cParams.method;
            obj.revolutions = cParams.revolutions;
            obj.timeSteps   = cParams.timeSteps;
        end

        function setHardCodedProperties(obj)
            x_lo=-1/2; x_up=1/2; y_lo=-1/2; y_up=1/2;
            obj.nElemX = 20;
            obj.nElemY = 20;
            [obj.coord,obj.connec] = CreateMesh(x_lo,x_up,y_lo,y_up,obj.nElemX,obj.nElemY);
            C      = zeros(size(obj.coord));
            C(:,1) = -obj.coord(:,2);
            C(:,2) = obj.coord(:,1);
            obj.a  = C;
        end
    end
end