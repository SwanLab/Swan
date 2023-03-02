% BPOPT Nonlinear Programming Solver
% Large-scale Interior Point solution methods
% Please note that this development version in MATLAB is primarily used
%   for algorithm development and small (<1000 variables) problems. An
%   optimized and large-scale version of this algorithm is available for
%   use at http://www.APMonitor.com and http://www.APOPT.com.
%
% Inputs
%   bp  = problem from bp_create
%
% Outputs
%   sol = solution returned as a structure
function [sol] = bpopt(bp)

% ensure that there is at least 1 input argument
if (nargin==0),
    disp('Insufficient arguments, call bp_create to create a problem');
    sol = [];
    return
end

disp('BPOPT Solver, MATLAB version')
disp('See http://apmonitor.com for more information')
disp('')
disp(['Solving Problem: ', bp.name])

% initial variables
[x,xL,xU,bL,bU] = bp_x_init(bp);
% add slack for inequality constraints only
k = 0;
m = max(size(bL));
si = [];
s = [];
sL = [];
sU = [];
for i = 1:m,
    if(bU(i)>bL(i)),
        k = k + 1;
        si(k) = i;
        s(k) = 0;
        sL(k) = bL(i);
        sU(k) = bU(i);
    end
end

% system size
% n variables
n = max(size(x));
% ns slack variables
ns = max(size(s));

% initial residuals
r = bp_res(bp,x,s,bL,bU);
% m constraints (inequality or equality)
m = max(size(r));
% ones vector
e = ones(n+ns,1);

% initial equation slack variables
k = 0;
for i = 1:m,
    % no slack variable when bU(i) == bL(i)
    if(bU(i)>bL(i)),
        k = k + 1;
        if (bp.slack_init),
           s(k) = r(i);  % corrected, was -r(i)
        else
           s(k) = 0.01;
        end
    end
end

% check for consistent bounds
if (min(xU-xL)<0),
    disp('bounds error (xU < xL)')
    sol = [];
    return
end
if (min(bU-bL)<0),
    disp('bounds error (bU < bL)')
    sol = [];
    return
end

% move x into feasible region and off boundary initially
for i = 1:n,
    if (x(i)<=xL(i) || x(i)>=xU(i)),
        fprintf(1,'Moving x0 to interior region\n')
        break
    end
end
% Move x variables to be feasible
for i = 1:n,
    if (x(i)<=xL(i)),
        x(i) = min(xU(i),xL(i)+1e-2);
    end
    if (x(i)>=xU(i)),
        x(i) = max(xL(i),xU(i)-1e-2);
    end
end
% Move slack variables to be feasible
for i = 1:ns,
    if (s(i)<=sL(i)),
        s(i) = min(sU(i),sL(i)+1e-2);
    end
    if (s(i)>=sU(i)),
        s(i) = max(sL(i),sU(i)-1e-2);
    end
end

tau = min(bp.tau_max,100*bp.mu);

% variable constraint multipliers
% zL*(x-xL) = mu  =>  zL = mu / (x-xL)
%zL(1:n) = 1;
for i = 1:n,
    zL(i) = bp.mu / (x(i)-xL(i));
end
for i = 1:ns,
    zL(n+i) = bp.mu / (s(i)-sL(i));
end
% zU*(xU-x) = mu  =>  zU = mu / (xU-x)
for i = 1:n,
    zU(i) = bp.mu / (xU(i)-x(i));
end
for i = 1:ns,
    zU(n+i) = bp.mu / (sU(i)-s(i));
end

% equation constraint multipliers initialization
g = bp_objgrad(bp,x,s);
J = bp_jac(bp,x,bL,bU);

lam = pinv(full(J*J'))*J*(zL'-zU'-g');

lam = lam';

% verify 1st and 2nd derivatives
if (bp.prob==0),
    ok_1st = true;
    ok_2nd = true;
else
    % test only for hard coded problems
    ok_1st = bp_verify_jac(bp,x,s,bL,bU);
    ok_2nd = bp_verify_hes(bp,x,s,lam,bL,bU);
    % test objective function gradients
    ok_obj_grad = bp_verify_objgrad(bp,x,s);
end

% initial parameters
alpha_pr = 1.0;
alpha_du = 1.0;

% 1-norm of infeasibilities
th = bp_theta(bp,x,s,bL,bU);
th_max = 10^4  * max(1,th);
th_min = 10^-4 * max(1,th);

% initialize iteration count
iter = 0;

% initialize filter
filter = [th_max bp_phi(bp,x,xL,xU,s,bL,bU)];

% print iteration zero
bp_iprint(bp,iter,x,lam,zL,zU,alpha_pr,alpha_du,s,xL,xU,bL,bU);

store = [iter x alpha_pr alpha_du bp.mu];

% start iterating
for iter = 1:bp.maxiter,
    
    % new residuals
    r = bp_res(bp,x,s,bL,bU);
    
    % 1-norm of infeasibilities
    th = bp_theta(bp,x,s,bL,bU);
    
    % phi, objective function and barrier terms
    ph = bp_phi(bp,x,xL,xU,s,bL,bU);
    
    if (bp.z_update==2),
        % update explicitly from z = mu / x
        for i = 1:n,
            zL(i) = bp.mu / (x(i)-xL(i));
        end
        for i = 1:ns,
            zL(n+i) = bp.mu / (s(i)-sL(i));
        end
        % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
        for i = 1:n,
            zU(i) = bp.mu / (xU(i)-x(i));
        end
        for i = 1:ns,
            zU(n+i) = bp.mu / (sU(i)-s(i));
        end
    end
    % make diagonal matrices
    zML = diag(zL);
    zMU = diag(zU);
    
    % sigmas
    dL = [x-xL s-sL];
    dU = [xU-x sU-s];
    dML = diag(dL);
    dMU = diag(dU);
    inv_dML = zeros(n+ns,n+ns);
    inv_dMU = zeros(n+ns,n+ns);
    SigL = zeros(n+ns,n+ns);
    SigU = zeros(n+ns,n+ns);
    for i = 1:n+ns,
       inv_dML(i,i) = 1 / dML(i,i); 
       inv_dMU(i,i) = 1 / dMU(i,i); 
       SigL(i,i) = zML(i,i) / dML(i,i);
       SigU(i,i) = zMU(i,i) / dMU(i,i);
    end
    
    % 1st and 2nd derivatives
    J = bp_jac(bp,x,bL,bU);
    W = bp_hes(bp,x,s,lam);
    g = bp_objgrad(bp,x,s);
    H = W + SigL + SigU;
    
    % Construct and solve linear system Ax=b
    
    % Symmetric, zL/zU explicitly solved
    % construct A and b
    A1 = [H,J';J,zeros(m,m)];
    %b1(1:n+ns,1) = g' - zL' + zU' + J'*lam';
    b1(1:n+ns,1) = g' - zL' + zU' + J'*lam' + zL'-bp.mu*inv_dML*e -zU'+bp.mu*inv_dMU*e;
    b1(n+ns+1:n+ns+m,1) = r(1:m)';
    
    % compute search direction, solving Ax = b
    d1 = -pinv(full(A1))*b1;
    
    dx1 = d1(1:n,1);
    ds1 = d1(n+1:n+ns,1);
    dlam1 = d1(n+ns+1:n+ns+m,1);
    % compute search direction for z (explicit solution)
    dzL1 = bp.mu*inv_dML*e - zL' - SigL*[dx1; ds1];
    dzU1 = bp.mu*inv_dMU*e - zU' + SigU*[dx1; ds1];
    search1 = [d1;dzL1;dzU1];

    % Un-symmetric, zL/zU implicitly solved
    % construct A
    A2 = [W J' -eye(n+ns) eye(n+ns);
        J zeros(m,m+2*(n+ns));
        diag(zL) zeros((n+ns),m) dML zeros(n+ns);
        -diag(zU) zeros(n+ns,m+n+ns) dMU];
    b2(1:n+ns,1) = g' - zL' + zU' + J'*lam';
    b2(n+ns+1:n+ns+m,1) = r(1:m)';
    b2(n+ns+m+1:2*(n+ns)+m) = dML*diag(zL)*e - bp.mu*e;
    b2(2*(n+ns)+m+1:3*(n+ns)+m) = dMU*diag(zU)*e - bp.mu*e;
    d2 = -pinv(full(A2))*b2;
    dx2 = d2(1:n,1);
    ds2 = d2(n+1:n+ns,1);
    dlam2 = d2(n+ns+1:n+ns+m,1);
    dzL2 = d2(n+ns+m+1:2*(n+ns)+m);
    dzU2 = d2(2*(n+ns)+m+1:3*(n+ns)+m);
    search2 = d2;

    % test for positive definiteness
    % A1 --------------------------------------------------
    % test for positive definiteness of A1
    if (bp.idebug>=2),
        if (min(eigs(A1))<0),
            fprintf(1,'A1 not pos definite - eigenvalue test\n');
        end
    end
    for i=1:length(A1),
        if ( det( A1(1:i, 1:i) ) <= 0 )
            isposdef = false;
            break;
        end
    end
    if(~isposdef && bp.idebug>=2),
        fprintf(1,'A1 not pos definite - determinant test\n');
    end    
    % A2 --------------------------------------------------
    % test for positive definiteness of A2
    isposdef = true;
    for i=1:length(A2),
        if ( det( A2(1:i, 1:i) ) <= 0 )
            isposdef = false;
            break;
        end
    end
    if(~isposdef && bp.idebug>=2),
        fprintf(1,'A2 not pos definite\n');
    end
    if (bp.idebug>=2),
        if (min(eigs(A2))<0),
            fprintf(1,'A2 not pos definite - eigenvalue test\n');
        end
    end
            
    if (bp.idebug>=2),
        fprintf(1,'Diff: %12.4e, Cond(A1): %12.4e, Cond(A2): %12.4e\n', ...
            sum(abs(search1-search2)),cond(full(A1)),cond(full(A2)));
    end
    
    % reduced or full matrix inversion
    %  1 = condensed, symmetric matrix
    %  2 = full, unsymmetric matrix
    if (bp.matrix==1),
        dx = dx1;
        ds = ds1;
        dlam = dlam1;
        dzL = dzL1;
        dzU = dzU1;
        if(cond(full(A1))>1e10 && bp.idebug>=2),
            fprintf(1,'Warning: A1 condition number high %12.4e\n',cond(full(A1)))
        end
    else
        dx = dx2;
        ds = ds2;
        dlam = dlam2;
        dzL = dzL2;
        dzU = dzU2;
        if(cond(full(A2))>1e10 && bp.idebug>=2),
            fprintf(1,'Warning: A2 condition number high %12.4e\n',cond(full(A2)))
        end
    end
        
    % compute acceptance point
    xa = x + alpha_pr * dx';
    if (ns>=1),
        sa = s + alpha_pr * ds';
    else
        sa = [];
    end
    lama = lam + alpha_du * dlam';
                
    % updating zLa and zUa
    switch(bp.z_update)
        case(1)
            % update from direct solve approach
            zLa = zL + alpha_du * dzL';
            zUa = zU + alpha_du * dzU';
        case(2)
            % update explicitly from z = mu / x
            for i = 1:n,
                zLa(i) = bp.mu / (xa(i)-xL(i));
                dzL(i,1) = zLa(i) - zL(i);
            end
            for i = 1:ns,
                zLa(n+i) = bp.mu / (sa(i)-sL(i));
                dzL(n+i,1) = zLa(n+i) - zL(n+i);
            end
            % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
            for i = 1:n,
                zUa(i) = bp.mu / (xU(i)-xa(i));
                dzU(i,1) = zUa(i) - zU(i);
            end
            for i = 1:ns,
                zUa(n+i) = bp.mu / (sU(i)-sa(i));
                dzU(n+i,1) = zUa(n+i) - zU(n+i);
            end
    end

    %xa
    %sa
    %lama
    %zLa
    %zUa
    
    % max alpha is that which brings the search point to within "tau" of constraint
    % tau is 0 to 0.01 (tau = mu when mu<0.01, otherwise tau=0.01)
    alpha_pr_max = 1.0;
    alpha_du_max = 1.0;
    % check for constraint violations
    for i = 1:n,
        if(xa(i)<xL(i)),
            alpha_pr_max = min(alpha_pr_max,(xL(i)+tau*(x(i)-xL(i))-x(i))/dx(i,1));
        end
        if(xa(i)>xU(i)),
            alpha_pr_max = min(alpha_pr_max,(xU(i)-tau*(xU(i)-x(i))-x(i))/dx(i,1));
        end
    end
    for i = 1:ns,
        if(sa(i)<sL(i)),
            alpha_pr_max = min(alpha_pr_max,(sL(i)+tau*(s(i)-sL(i))-s(i))/ds(i,1));
        end
        if(sa(i)>sU(i)),
            alpha_pr_max = min(alpha_pr_max,(sU(i)+tau*(sU(i)-s(i))-s(i))/ds(i,1));
        end
    end
    for i = 1:n,
        if (bp.z_update==1),            
            if(zLa(i)<0),
                alpha_du_max = min(alpha_du_max,(tau*zL(i)-zL(i))/dzL(i,1));
                %zLa(i) = 0;
            end
            if(zUa(i)<0),
                alpha_du_max = min(alpha_du_max,(tau*zU(i)-zU(i))/dzU(i,1));
                %zUa(i) = 0;
            end
        end
    end
    
    
    % line search
    %  1 = reduction in merit function
    %  2 = simple clipping
    %  3 = filter method
    switch bp.line_search,
        case(1)
            % alpha_du set as approach to constraint
            alpha_du = alpha_du_max;
            % compute new acceptance point
            alpha_pr = min(alpha_pr,alpha_pr_max);
        case(2)
            % test case for alpha_pr, alpha_du = 1 with clipping only
            alpha_pr = 1;
            alpha_du = 1;
        case(3)
            % alpha_du set as approach to constraint
            alpha_du = alpha_du_max;
            % compute new acceptance point
            alpha_pr = min(alpha_pr,alpha_pr_max);
    end
    
    % apply alpha values
    xa = x + alpha_pr * dx';
    if(ns>=1),
        sa = s + alpha_pr * ds';
    end
    lama = lam + alpha_pr * dlam';

    % updating zLa and zUa
    switch(bp.z_update)
        case(1)
            % update from direct solve approach
            zLa = zL + alpha_du * dzL';
            zUa = zU + alpha_du * dzU';
        case(2)
            % update explicitly from z = mu / x
            for i = 1:n,
                zLa(i) = bp.mu / (xa(i)-xL(i));
                dzL(i,1) = zLa(i) - zL(i);
            end
            for i = 1:ns,
                zLa(n+i) = bp.mu / (sa(i)-sL(i));
                dzL(n+i,1) = zLa(n+i) - zL(n+i);
            end
            % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
            for i = 1:n,
                zUa(i) = bp.mu / (xU(i)-xa(i));
                dzU(i,1) = zUa(i) - zU(i);
            end
            for i = 1:ns,
                zUa(n+i) = bp.mu / (sU(i)-sa(i));
                dzU(n+i,1) = zUa(n+i) - zU(n+i);
            end
    end
    
    % clipping (this should already be arranged by alpha_max)
    % push away from the boundary with tau
    for i = 1:n,
        if(xa(i)<xL(i)),
            xa(i) = xL(i)+tau*(x(i)-xL(i));
        end
        if(xa(i)>xU(i)),
            xa(i) = xU(i)-tau*(xU(i)-x(i));
        end
        if(zLa(i)<0),
            zLa(i) = tau*(zL(i));
        end
        if(zUa(i)<0),
            zUa(i) = tau*(zU(i));
        end
    end
    for i = 1:ns,
        if(sa(i)<sL(i)),
            sa(i) = sL(i)+tau*(s(i)-sL(i));
        end
        if(sa(i)>sU(i)),
            sa(i) = sU(i)-tau*(sU(i)-s(i));
        end
    end
    
    % predicted reduction in the merit function
    pred = -alpha_pr*g*[dx;ds] - 0.5*alpha_pr^2*[dx;ds]'*H*[dx;ds] + bp.nu*(norm(r',1)-norm(r'+alpha_pr*J*[dx;ds],1));
    ared = bp_merit(bp,x,xL,xU,s,bL,bU) - bp_merit(bp,xa,xL,xU,sa,bL,bU);

    % line search criteria
    %  1 = reduction in merit function
    %  2 = simple clipping
    %  3 = filter method
    switch bp.line_search,
        case(1)
            % merit function
            pred_phi_1 = g*[dx;ds];
            % set 2nd derivative contribution to zero if < 0
            pred_phi_2 = max(0,0.5*[dx;ds]'*H*[dx;ds]);
            pred_phi_decrease = pred_phi_1 + pred_phi_2;
            th_a = bp_theta(bp,xa,sa,bL,bU);
            rho = 0.1;
            new_nu = pred_phi_decrease / ((1-rho) * th_a);
            %if (new_nu>bp.nu),
            %    bp.nu = min(1000,new_nu + 1);
            %end
            bp.nu = max(1,min(1000,new_nu));
            % update predicted and actual reductions with new nu value
            pred = -alpha_pr*g*[dx;ds] - 0.5*alpha_pr^2*[dx;ds]'*H*[dx;ds] + bp.nu*(norm(r',1)-norm(r'+alpha_pr*J*[dx;ds],1));
            ared = bp_merit(bp,x,xL,xU,s,bL,bU) - bp_merit(bp,xa,xL,xU,sa,bL,bU);
            eta = 0.2;
            % compare actual reduction to predicted reduction
            % as long as the actual reduction is a fraction of the
            % predicted reduction then accept the trial point
            if (ared>=eta*pred),
                ac = true;
            else
                ac = false;
            end
        case(2)
            ac = true;
            % test case for alpha_pr, alpha_du = 1 with clipping only
            alpha_pr = 1;
            alpha_du = 1;
        case(3)
            % check if acceptable point with filter method
            [ac] = bp_accept(bp,x,xa,xL,xU,filter);
    end
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing for acceptance
    %ac = true;
    % second order correction
    if (false),
        if (not(ac)),
            % new residuals
            b1(n+1:n+m,1) = alpha_pr*bp_res(bp,x,s,bL,bU)' + bp_res(bp,xa,sa,bL,bU)';
            % compute search direction, solving Ax = b
            d3 = -pinv(full(A1))*b1;
            dx3 = d3(1:n,1);
            ds3 = d3(n+1:n+ns,1);
            dlam3 = d3(n+ns+1:n+ns+m,1);
            % compute search direction for z (explicit solution)
            dzL3 = bp.mu*inv_dML*e - zL' - SigL*[dx3; ds3];
            dzU3 = bp.mu*inv_dMU*e - zU' + SigU*[dx3; ds3];
            
            % compute new acceptance point
            xa = x + alpha_pr * dx3';
            sa = s + alpha_pr * ds3';
            lama = lam + alpha_du * dlam3';
            zLa = zL + alpha_du * dzL3';
            zUa = zU + alpha_du * dzU3';
            
            % check for acceptance of trial point
            [ac] = bp_accept(bp,x,xa,xL,xU,filter);
            if (ac),
                fprintf(1,'2nd order correction success\n');
            else
                fprintf(1,'2nd order correction failed to find accepatable point\n');
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % optionally show contour plot
    if (bp.contour),
       bp_contour(bp,x,alpha_pr * dx);
    end
    
    if (ac),
        % accept point
        x = xa;
        s = sa;
        lam = lama;
        zL = zLa;
        zU = zUa;
        
        % update filter
        filter = [filter; th ph];
        %else
        %if (bp_theta(bp,xa)>bp_theta(bp,x)),
        % apply second order correction
        
        %end
    end
    
    % check for convergence
    s_max = 100; % > 1
    s_d = max(s_max,(sum(abs(lam))+sum(abs(zL))+sum(abs(zU)))/(m+2*(n+ns)));
    s_c = max(s_max,(sum(abs(zU))+sum(abs(zL)))/(2*(n+ns)));
    
    part(1) = max(abs(bp_objgrad(bp,x,s)' + bp_jac(bp,x,bL,bU)'*lam' - zL' + zU'))/s_d;
    part(2) = max(abs(bp_res(bp,x,s,bL,bU)));
    % use mu = 0 to test for convergence
    %part(3) = max(abs(diag([x-xL s-sL])*diag(zL)*e - bp.mu*e))/s_c;
    %part(4) = max(abs(diag([xU-x sU-s])*diag(zU)*e - bp.mu*e))/s_c;
    part(3) = max(abs(diag([x-xL s-sL])*diag(zL)*e))/s_c;
    part(4) = max(abs(diag([xU-x sU-s])*diag(zU)*e))/s_c;
    e_mu = max(part);
    
    if (bp.idebug>=1),
        fprintf(1,'pred: %12.4e, ared: %12.4e, e_mu: %12.4e, k*mu: %12.4e\n',pred,ared,e_mu,0.2*bp.mu);
    end
    
    % check for termination conditions
    if (e_mu <= bp.e_tol),
        fprintf(1,'\nSuccessful solution\n');
        status = 'success';
        break;
    end
    
    % check for new barrier problem
    k_mu  = 0.2; % (0,1)
    
    % only update mu if requested (not debugging)
    if (bp.mu_update),
        if (e_mu < k_mu * bp.mu),
            th_mu = 1.5; % (1,2)
            % update mu
            bp.mu = max(bp.e_tol/10,min(k_mu*bp.mu,bp.mu^th_mu));
            % update tau
            tau = min(bp.tau_max,100*bp.mu);
            % re-initialize filter
            filter = [th_max bp_phi(bp,x,xL,xU,s,bL,bU)];
        end
    end
    % print iteration
    bp_iprint(bp,iter,x,lam,zL,zU,alpha_pr,alpha_du,s,xL,xU,bL,bU);
    
    store = [store; iter x alpha_pr alpha_du bp.mu];
    
    % reached the end of the iteration loop without convergence
    if (iter==bp.maxiter),
        status = 'failed: max iterations';
        break
    end
    
    % don't do a line search in this MATLAB version
    % just cycle through on another iteration with a lower alpha if
    %   the point was not accepted
    if (ac),
        % reset alpha
        alpha_pr = 1.0;
        alpha_du = 1.0;
    else
        % reject point and move alpha_x
        alpha_pr = alpha_pr / 2;
        alpha_du = alpha_du / 2;
        if (alpha_pr < 1e-4),
            alpha_pr = 1.0;
            alpha_du = 1.0;
        end
    end
    
end

%% Display final solution
disp('Status')
disp(status)
disp('Solution: ')
disp(x)
disp('Slack Variables: ')
disp(s)
disp('Equation Multipliers: ')
disp(lam)
disp('Lower Constraint Variable Mult: ')
disp(zL)
disp('Upper Constraint Variable Mult: ')
disp(zU)

if (bp.prob>=1),
    %% Display figures on iteration progress
    figure(1)
    hold off;
    i = 1;
    plot(store(:,1),store(:,i+1),'k-')
    hold on;
    i = 2;
    plot(store(:,1),store(:,i+1),'b-')
    legend('x_1','x_2')
    
    figure(2)
    hold off;
    plot(store(:,1),store(:,n+2),'r-')
    hold on;
    plot(store(:,1),store(:,n+3),'b-')
    plot(store(:,1),log10(store(:,n+4)),'g-')
    legend('alpha_{pr}','alpha_{du}','log_{10}(mu)');
end

% record solution for function return
sol.status = status;
sol.x = x;
sol.s = s;
sol.lam = lam;
sol.zL = zL;
sol.zU = zU;
