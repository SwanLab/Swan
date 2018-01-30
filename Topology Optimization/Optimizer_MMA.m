classdef Optimizer_MMA < Optimizer
    properties
        kkttol
        kktnorm
        maxoutit
        x
        xold1
        xold2
        xmin
        xmax
        low
        outit=0;
        outeriter=0;
        upp
        m
        c
        d
        a0
        a
        n
        f0val
        df0dx
        fval
        dfdx
    end 
    methods
        function obj=Optimizer_MMA(settings)
            obj@Optimizer(settings);           
            obj.maxoutit=1e4;
        end
        function kkttol=get.kkttol(obj)
            kkttol=obj.target_parameters.optimality_tol;
        end          
        function x=updateX(obj,x,cost,constraint,interpolation,filter)      
                obj.checkInitial(x,cost,constraint);
                obj.outit=obj.outit+1;
                iter=obj.outit
                obj.outeriter = obj.outeriter+1;
                %%%% The MMA subproblem is solved at the point xval:
                [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,obj.low,obj.upp] = ...
                    obj.mmasub(obj.m,obj.n,obj.outeriter,x,obj.xmin,obj.xmax,obj.xold1,obj.xold2, ...
                    obj.f0val,obj.df0dx,obj.fval,obj.dfdx,obj.low,obj.upp,obj.a0,obj.a,obj.c,obj.d);
                %%%% Some vectors are updated:
                obj.xold2 = obj.xold1;
                obj.xold1 = x;
                x = xmma;
                %%%% The user should now calculate function values and gradients
                %%%% of the objective- and constraint functions at xval.
                %%%% The results should be put in f0val, df0dx, fval and dfdx.
                obj.update_physical_variables(x,interpolation,filter); 
                cost.computef(x,obj.physicalProblem,interpolation,filter);
                constraint.computef(x,obj.physicalProblem,interpolation,filter);

                [obj.f0val,obj.df0dx,obj.fval,obj.dfdx] = obj.funmma(cost,constraint);
                %%%% The residual vector of the KKT conditions is calculated:
                [~,obj.kktnorm] = obj.kktcheck(obj.m,obj.n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
                    obj.xmin,obj.xmax,obj.df0dx,obj.fval,obj.dfdx,obj.a0,obj.a,obj.c,obj.d);
   
                obj.stop_criteria = obj.kktnorm > obj.kkttol && obj.outit < obj.maxoutit;
        end
        function checkInitial(obj,x_ini,cost,constraint)
            if isempty(obj.x)
            obj.x = x_ini;
            obj.xold1=obj.x;
            obj.xold2 = obj.xold1;
            obj.xmin=zeros(length(x_ini),1);
            obj.xmax=ones(length(x_ini),1);
            obj.low = obj.xmin;
            obj.upp = obj.xmax;           
            [obj.f0val,obj.df0dx,obj.fval,obj.dfdx]=obj.funmma(cost, constraint);
            obj.m = length(obj.fval);
            obj.c= 1e3*ones(obj.m,1);
            obj.d = 0*ones(obj.m,1);
            obj.a0 = 1;
            obj.a = 0*ones(obj.m,1);
            obj.n = length(obj.x);
            end
        end
        function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
                mmasub(obj,m,n,iter,xval,xmin,xmax,xold1,xold2, ...
                ~,df0dx,fval,dfdx,low,upp,a0,a,c,d)
            epsimin = 10^(-7);
            raa0 = 1e-5;
            move = 1.0;
            albefa = 0.1;
            asyinit = 0.5;
            asyincr = 1.2;
            asydecr = 0.7;
            eeen = ones(n,1);
            eeem = ones(m,1);
            zeron = zeros(n,1);
            
            % Calculation of the asymptotes low and upp :
            if iter < 2.5
                low = xval - asyinit*(xmax-xmin);
                upp = xval + asyinit*(xmax-xmin);
            else
                zzz = (xval-xold1).*(xold1-xold2);
                factor = eeen;
                factor(find(zzz > 0)) = asyincr;
                factor(find(zzz < 0)) = asydecr;
                low = xval - factor.*(xold1 - low);
                upp = xval + factor.*(upp - xold1);
                lowmin = xval - 10*(xmax-xmin);
                lowmax = xval - 0.01*(xmax-xmin);
                uppmin = xval + 0.01*(xmax-xmin);
                uppmax = xval + 10*(xmax-xmin);
                low = max(low,lowmin);
                low = min(low,lowmax);
                upp = min(upp,uppmax);
                upp = max(upp,uppmin);
            end
            
            % Calculation of the bounds alfa and beta :
            
            zzz1 = low + albefa*(xval-low);
            zzz2 = xval - move*(xmax-xmin);
            zzz  = max(zzz1,zzz2);
            alfa = max(zzz,xmin);
            zzz1 = upp - albefa*(upp-xval);
            zzz2 = xval + move*(xmax-xmin);
            zzz  = min(zzz1,zzz2);
            beta = min(zzz,xmax);
            
            % Calculations of p0, q0, P, Q and b.
            
            xmami = xmax-xmin;
            xmamieps = 1e-5*eeen;
            xmami = max(xmami,xmamieps);
            xmamiinv = eeen./xmami;
            ux1 = upp-xval;
            ux2 = ux1.*ux1;
            xl1 = xval-low;
            xl2 = xl1.*xl1;
            uxinv = eeen./ux1;
            xlinv = eeen./xl1;
            %
            p0 = zeron;
            q0 = zeron;
            p0 = max(df0dx,0);
            q0 = max(-df0dx,0);
            %p0(find(df0dx > 0)) = df0dx(find(df0dx > 0));
            %q0(find(df0dx < 0)) = -df0dx(find(df0dx < 0));
            pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
            p0 = p0 + pq0;
            q0 = q0 + pq0;
            p0 = p0.*ux2;
            q0 = q0.*xl2;
            %
            dfdx = sparse(dfdx);
            P = max(dfdx,0);
            Q = max(-dfdx,0);
            %P(find(dfdx > 0)) = dfdx(find(dfdx > 0));
            %Q(find(dfdx < 0)) = -dfdx(find(dfdx < 0));
            PQ = 0.001*(P + Q) + raa0*eeem*xmamiinv';
            %PQ = 0.001*(P + Q) + raa0*eeem*xmamiinv;
            P = P + PQ;
            Q = Q + PQ;
            P = P * spdiags(ux2,0,n,n);
            Q = Q * spdiags(xl2,0,n,n);
            b = P*uxinv + Q*xlinv - fval ;
            %
            %%% Solving the subproblem by a primal-dual Newton method
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
                obj.subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
            
        end
    end
    methods (Static)
        function [f,df,c,dc] = funmma(cost,constraint)
            f=cost.value;
            df=cost.gradient;
            c=constraint.value;
            dc=constraint.gradient;
            dc = dc';
            
            %% Check constraint case
            %             if strcmp(constraint_type,'equality')
            %                 c = [c;-c];
            %                 dc = [dc;-dc];
            %             end
            
            %% Re-scale constraints
            % In many applications, the constraints are on the form yi(x) <= ymaxi
            % The user should then preferably scale the constraints in such a way that 1 <= ymaxi <= 100 for each i
            % (and not ymaxi = 10^10 for example).
            kconstr = 1;
            cconstr = 0;
            c = kconstr*c;
            c(c > 0) = c(c > 0) + cconstr;
            c(c < 0) = c(c < 0) - cconstr;
            % dc = kconstr*dc;
            
            %% Re-scale objective function
            % The objective function f(x) should preferably be scaled such that
            % 1 <= f0(x) <= 100 for reasonable values on the variables.
            kfun = 1;
            cfun = 0;
            f = kfun*f + cfun;
            df = kfun*df;
            
        end
        function [xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma] = ...
                subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
            %
            % This function subsolv solves the MMA subproblem:
            %
            % minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
            %          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
            %
            % subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
            %            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
            %
            % Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
            % Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
            %
            een = ones(n,1);
            eem = ones(m,1);
            epsi = 1;
            epsvecn = epsi*een;
            epsvecm = epsi*eem;
            x = 0.5*(alfa+beta);
            y = eem;
            z = 1;
            lam = eem;
            xsi = een./(x-alfa);
            xsi = max(xsi,een);
            eta = een./(beta-x);
            eta = max(eta,een);
            mu  = max(eem,0.5*c);
            zet = 1;
            s = eem;
            itera = 0;
            while epsi > epsimin
                epsvecn = epsi*een;
                epsvecm = epsi*eem;
                ux1 = upp-x;
                xl1 = x-low;
                ux2 = ux1.*ux1;
                xl2 = xl1.*xl1;
                uxinv1 = een./ux1;
                xlinv1 = een./xl1;
                plam = p0 + P'*lam ;
                qlam = q0 + Q'*lam ;
                gvec = P*uxinv1 + Q*xlinv1;
                dpsidx = plam./ux2 - qlam./xl2 ;
                rex = dpsidx - xsi + eta;
                rey = c + d.*y - mu - lam;
                rez = a0 - zet - a'*lam;
                relam = gvec - a*z - y + s - b;
                rexsi = xsi.*(x-alfa) - epsvecn;
                reeta = eta.*(beta-x) - epsvecn;
                remu = mu.*y - epsvecm;
                rezet = zet*z - epsi;
                res = lam.*s - epsvecm;
                residu1 = [rex' rey' rez]';
                residu2 = [relam' rexsi' reeta' remu' rezet res']';
                residu = [residu1' residu2']';
                residunorm = sqrt(residu'*residu);
                residumax = max(abs(residu));
                ittt = 0;
                while residumax > 0.9*epsi & ittt < 200
                    ittt=ittt + 1;
                    itera=itera + 1;
                    ux1 = upp-x;
                    xl1 = x-low;
                    ux2 = ux1.*ux1;
                    xl2 = xl1.*xl1;
                    ux3 = ux1.*ux2;
                    xl3 = xl1.*xl2;
                    uxinv1 = een./ux1;
                    xlinv1 = een./xl1;
                    uxinv2 = een./ux2;
                    xlinv2 = een./xl2;
                    plam = p0 + P'*lam ;
                    qlam = q0 + Q'*lam ;
                    gvec = P*uxinv1 + Q*xlinv1;
                    GG = P*spdiags(uxinv2,0,n,n) - Q*spdiags(xlinv2,0,n,n);
                    dpsidx = plam./ux2 - qlam./xl2 ;
                    delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x);
                    dely = c + d.*y - lam - epsvecm./y;
                    delz = a0 - a'*lam - epsi/z;
                    dellam = gvec - a*z - y - b + epsvecm./lam;
                    diagx = plam./ux3 + qlam./xl3;
                    diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x);
                    diagxinv = een./diagx;
                    diagy = d + mu./y;
                    diagyinv = eem./diagy;
                    diaglam = s./lam;
                    diaglamyi = diaglam+diagyinv;
                    if m < n
                        blam = dellam + dely./diagy - GG*(delx./diagx);
                        bb = [blam' delz]';
                        Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
                        AA = [Alam     a
                            a'    -zet/z ];
                        solut = AA\bb;
                        dlam = solut(1:m);
                        dz = solut(m+1);
                        dx = -delx./diagx - (GG'*dlam)./diagx;
                    else
                        diaglamyiinv = eem./diaglamyi;
                        dellamyi = dellam + dely./diagy;
                        Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
                        azz = zet/z + a'*(a./diaglamyi);
                        axz = -GG'*(a./diaglamyi);
                        bx = delx + GG'*(dellamyi./diaglamyi);
                        bz  = delz - a'*(dellamyi./diaglamyi);
                        AA = [Axx   axz
                            axz'  azz ];
                        bb = [-bx' -bz]';
                        solut = AA\bb;
                        dx  = solut(1:n);
                        dz = solut(n+1);
                        dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
                    end
                    %
                    dy = -dely./diagy + dlam./diagy;
                    dxsi = -xsi + epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa);
                    deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x);
                    dmu  = -mu + epsvecm./y - (mu.*dy)./y;
                    dzet = -zet + epsi/z - zet*dz/z;
                    ds   = -s + epsvecm./lam - (s.*dlam)./lam;
                    xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']';
                    dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
                    %
                    stepxx = -1.01*dxx./xx;
                    stmxx  = max(stepxx);
                    stepalfa = -1.01*dx./(x-alfa);
                    stmalfa = max(stepalfa);
                    stepbeta = 1.01*dx./(beta-x);
                    stmbeta = max(stepbeta);
                    stmalbe  = max(stmalfa,stmbeta);
                    stmalbexx = max(stmalbe,stmxx);
                    stminv = max(stmalbexx,1);
                    steg = 1/stminv;
                    %
                    xold   =   x;
                    yold   =   y;
                    zold   =   z;
                    lamold =  lam;
                    xsiold =  xsi;
                    etaold =  eta;
                    muold  =  mu;
                    zetold =  zet;
                    sold   =   s;
                    %
                    itto = 0;
                    resinew = 2*residunorm;
                    while resinew > residunorm & itto < 50
                        itto = itto+1;
                        x   =   xold + steg*dx;
                        y   =   yold + steg*dy;
                        z   =   zold + steg*dz;
                        lam = lamold + steg*dlam;
                        xsi = xsiold + steg*dxsi;
                        eta = etaold + steg*deta;
                        mu  = muold  + steg*dmu;
                        zet = zetold + steg*dzet;
                        s   =   sold + steg*ds;
                        ux1 = upp-x;
                        xl1 = x-low;
                        ux2 = ux1.*ux1;
                        xl2 = xl1.*xl1;
                        uxinv1 = een./ux1;
                        xlinv1 = een./xl1;
                        plam = p0 + P'*lam ;
                        qlam = q0 + Q'*lam ;
                        gvec = P*uxinv1 + Q*xlinv1;
                        dpsidx = plam./ux2 - qlam./xl2 ;
                        rex = dpsidx - xsi + eta;
                        rey = c + d.*y - mu - lam;
                        rez = a0 - zet - a'*lam;
                        relam = gvec - a*z - y + s - b;
                        rexsi = xsi.*(x-alfa) - epsvecn;
                        reeta = eta.*(beta-x) - epsvecn;
                        remu = mu.*y - epsvecm;
                        rezet = zet*z - epsi;
                        res = lam.*s - epsvecm;
                        residu1 = [rex' rey' rez]';
                        residu2 = [relam' rexsi' reeta' remu' rezet res']';
                        residu = [residu1' residu2']';
                        resinew = sqrt(residu'*residu);
                        steg = steg/2;
                    end
                    residunorm=resinew;
                    residumax = max(abs(residu));
                    steg = 2*steg;
                end
                if ittt > 198
                    epsi
                    ittt
                end
                epsi = 0.1*epsi;
            end
            xmma   =   x;
            ymma   =   y;
            zmma   =   z;
            lamma =  lam;
            xsimma =  xsi;
            etamma =  eta;
            mumma  =  mu;
            zetmma =  zet;
            smma   =   s;
        end
        function [residu,residunorm,residumax] = ...
                kktcheck(~,~,x,y,z,lam,xsi,eta,mu,zet,s, ...
                xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
            
            %
            rex   = df0dx + dfdx'*lam - xsi + eta;
            rey   = c + d.*y - mu - lam;
            rez   = a0 - zet - a'*lam;
            relam = fval - a*z - y + s;
            rexsi = xsi.*(x-xmin);
            reeta = eta.*(xmax-x);
            remu  = mu.*y;
            rezet = zet*z;
            res   = lam.*s;
            %
            residu1 = [rex' rey' rez]';
            residu2 = [relam' rexsi' reeta' remu' rezet res']';
            residu = [residu1' residu2']';
            residunorm = sqrt(residu'*residu);
            residumax = max(abs(residu));
            %----------------------------------------
        end
    end
end