classdef RadialQuadrature
    % Approximation of a function as a series of the functions (e_i)
    % where e_i is the i-th normalized (in H10 norm) eigenfunction of the
    % Laplace operator with Dirichlet boundary conditions.
    % The approximation takes place on the set a(1) < |x| < a(2) of R^2
    % and writes
    % func \approx \sum_{i \in I} alpha0(i) e_i(x)
    % The coefficients $\alpha$ are chosen as the minimizers of the H10
    % norm of this approximation.
    properties (GetAccess = public, SetAccess = protected)
        a, b, kernel, tol, Pmax, alpha0, alpha, rho,
        resH10, errLinfStored, reachedTol, times, scal01, nIter;
        
    end
    
    methods
        %% constructor and displays
        
        function[rq] = RadialQuadrature(aa,kernel,ttol,varargin)
            if nargin > 0
                %% Input description
                % REQUIRED INPUTS :
                % - aa : as described in object description
                % - func : the radial function to be approximated
                % - der : the derivative of func
                % OPTIONAL INPUTS (Name-Value pairs)
                % - Pmax / type : integer or Inf / Default : Inf
                % Descr : maximal number of elements in quadrature
                % - verbose / type : logical / default : false
                % Descr : Set to true to enable text displays during computation
                % - batch_size : sets the constant K such that fix(K/a)+1 is the
                % number of frequencies added at each iteration.
                
                %% Creation of quadrature
                
                out = radialQuad(aa,kernel,ttol,varargin{:});
                
                %% Export
                rq.a = out.a;
                rq.b = out.b;
                rq.kernel = kernel;
                rq.tol = out.tol;
                rq.Pmax = out.Pmax;
                rq.alpha0 = out.alpha0;
                rq.alpha = out.alpha;
                rq.rho = out.rho;
                rq.resH10 = out.resH10;
                rq.errLinfStored = out.errLinf;
                rq.reachedTol = out.reachedTol;
                rq.times = out.times;
                rq.scal01 = out.scal01;
                rq.nIter = out.nIter;
                rrho = rq.rho;
                assert(rrho(1) == 0);
                % This convention must always be satisfied !! The constant.
            else
                % Represents the null function
                rq.a = 0;
                rq.b = Inf;
                rq.kernel = Kernel;
                rq.tol = 0;
                rq.Pmax = NaN;
                rq.alpha0 = 0;
                rq.alpha = 0;
                rq.rho = 0;
                rq.errLinfStored = [];
                rq.reachedTol = true;
                times = struct;
                times.freq = 0;
                times.storeValsError = 0;
                times.chol = 0;
                times.cholInv = 0;
                times.projections = 0;
                times.testErr = 0;
                times.total = 0;
                rq.times = times;
                rq.nIter = 0;
            end
        end
        function[] = show(this,n)
            noDer = false;
            if nargin==1
                showAux(this,'all',noDer);
            else
                showAux(this,n,noDer);
            end
        end
        function[] = showDer(this,n)
            optDer = true;
            if nargin==1
                showAux(this,'all',optDer);
            else
                showAux(this,n,optDer);
            end
        end
        function[] = disp(this)
            dispAux(this);
        end
        function[] = showTimes(this)
            t = this.times;
            fprintf('Time history (s) :\n\n')
            rowNames = {'J0 roots';'J0 values';'Cholesky';'Cholesky^{-1}';'Projection';'Linf';'TOTAL'};
            Times = struct2array(t)';
            T = table(Times,'RowNames',rowNames,'VariableNames',{'t'});
            disp(T);
            
            if sum(Times(1:end-1))>0
                explode = Times(1:end-1)*0+1;
                figure('Name','Computation time','NumberTitle','off')
                pie(Times(1:end-1)/sum(Times(1:end-1)),explode,rowNames(1:end-1))
            end
        end
        
        %% Other methods
        function[err] = errLinf(this)
            if isempty(this.errLinfStored)
                t = RadialQuadrature.tTestLinf(this.a,this.b,this.rho);
                [~,errs] = this.eval(t);
                err = max(abs(errs));
            else
                err = this.errLinfStored;
            end
        end
        function[s] = getQuadErrString(this)
            if this.errLinf < 0.1*this.tol
                s = sprintf('Quadrature error : < %s \n',num2str(this.errLinf));
            elseif this.errLinf <= this.tol
                s= sprintf('Quadrature error : <= %s \n',num2str(this.tol));
            else
                s = sprintf('Quadrature error : %s \n',num2str(this.errLinf));
                s = [s,sprintf('The tolerance defined by user was not reached. Increase Pmax \n')];
            end
        end
        function[rq] = plus(rq1,rq2)
            rq = RadialQuadrature;
            rq.a = max(rq1.a,rq2.a);
            rq.b = min(rq1.b,rq2.b);
            rq.kernel = rq1.kernel + rq2.kernel;
            rq.tol = rq1.tol + rq2.tol;
            rq.Pmax = max(rq1.Pmax,rq2.Pmax);
            [rq.rho,rq.alpha] = mergeQuads(rq1.rho,rq1.alpha,rq2.rho,rq2.alpha);
            [~,rq.alpha0] = mergeQuads(rq1.rho,rq1.alpha0,rq2.rho,rq2.alpha0);
            [~,rq.scal01] = mergeQuads(rq1.rho(2:end),rq1.scal01,rq2.rho(2:end),rq2.scal01);
            rq.resH10 = rq1.resH10 + rq2.resH10;
            rq.errLinfStored = rq.errLinf;
            rq.reachedTol = and(rq1.reachedTol,rq2.reachedTol);
            rq.times = sumStruct(rq1.times,rq2.times);
            rq.nIter = max(rq1.nIter,rq2.nIter);
        end
        function[rq] = mtimes(rq1,lambda)
            if isa(rq1,'RadialQuadrature')
                assert(and(isa(lambda,'double'),isscalar(lambda)));
                rq = rq1;
                rq.alpha = lambda*rq.alpha;
                rq.alpha0 = lambda*rq.alpha0;
                rq.kernel = lambda*rq.kernel;
                rq.tol = abs(lambda)*rq.tol;
                rq.resH10 = abs(lambda)*rq.resH10;
                rq.errLinfStored = abs(lambda)*rq1.errLinf;
                rq.scal01 = lambda*rq.scal01;
            else
                rq = mtimes(lambda,rq1);
            end
        end
        function[res] = getTheoreticalBoundH10(this)
            res =  sqrt(-log(this.a)/(2*pi))*this.resH10;
        end
        function[s] = getFunc(this)
            s = func2str(this.kernel.func);
        end
        function[this] = dilatation(this,lambda)
            assert(lambda > this.a); % Otherwise new value of a would be > 1
            this.a = this.a/lambda;
            this.b = min(this.b/lambda,1);
            oldRho = this.rho;
            this.rho = this.rho*lambda;
            newKernel = this.kernel.dilatation(lambda);
            this.kernel = newKernel;
            this.alpha0 = this.alpha0.*Cp(oldRho)./Cp(this.rho);
        end 
        function[quad,err] = eval(this,r)
            quad = coeff2func(this.alpha0,this.rho,r).';
            trueVal = this.kernel.func(r);
            err = quad - trueVal;
        end
        function[quad,err] = evalDer(this,r)
            quad = coeff2der(this.alpha0,this.rho,r);
            quad = quad(:);
            if nargout >=2
                trueVal = this.kernel.der(r);
                trueVal = trueVal(:);
                err = quad - trueVal;
            end
        end
    end
    methods (Static)
        function[t] = tTestLinf(a,b,rho)
            if a==0
                t = [linspace(0,0.1,fix(max(rho))+1),linspace(0.9,1,fix(max(rho))+1)];
            else
                t = [linspace(a,min(2*a,b),fix(10*max(rho)*min(a,b-a))+1)...
                    linspace(min(b,1-a),1,fix(10*max(rho)*min(1-b,a))+1)];  % 10 points per freq
            end
        end
    end
end


%% Auxiliary functions

function[] = dispAux(in)

func = in.kernel.func;
rho = in.rho;
a = in.a;
b = in.b;
P = length(rho);

if isa(in,'RadialQuadratureY0')
    s = 'Y0(Rx)';
    fprintf('Radial Quadrature of function %s \n',s);
    fprintf('With R = %s \n',num2str(in.getR));
elseif isa(in,'RadialQuadratureLog')
    s = 'log(Rx)';
    fprintf('Radial Quadrature of function %s \n',s);
    fprintf('With R = %s \n',num2str(in.getR));
else
    fprintf('Radial Quadrature of function %s \n',func2str(func));
end
fprintf('Domain of approximation : [%s, %s]\n',num2str(a),num2str(b))
fprintf('Number of components (not including constant) : %d \n',P-1)
disp(in.getQuadErrString);

end
function [] = showAux(in,n,optDer)
if nargin < 2
    n = 'all';
    optDer = false;
elseif nargin==2
    optDer = false;
end
a = in.a;
b = in.b;
der = in.kernel.der;
func = in.kernel.func;
rho = in.rho;
alpha0 = in.alpha0;
tol = in.tol;
scal01 = in.scal01;
startFreq = in.kernel.startFreq;

crop = 0.999; % Needed to avoid having a 10^-16 value in the error at the outer edge
% Since function and approx are always equal at b



if isequal(n,'all')
    if optDer
        figure('Name','Derivative of Radial quadrature','NumberTitle','off')
    else
        figure('Name','Radial quadrature','NumberTitle','off')
    end
    
    subplot(1,3,1)
    showAux(in,1,optDer);
    
    subplot(1,3,2)
    showAux(in,2,optDer);
    
    subplot(1,3,3)
    showAux(in,3,optDer)
else
    switch n
        case 1
            t2 = linspace(0,1,fix(min(max(rho)*5,20000)));
            t1 = t2(and(t2>a/3,t2<b*crop));
            if optDer 
                quad = in.evalDer(t2);
                plot(t1,real(der(t1)),'b');
                hold on
                plot(t1,imag(der(t1)),'r');
            else
                quad = in.eval(t2);
                plot(t1,real(func(t1)),'b');
                hold on
                plot(t1,imag(func(t1)),'r');
            end
            plot(t2,real(quad),'b--');
            plot(t2,imag(quad),'r--');
            plot([a a],ylim,'k--');
            if b<1
                plot([a a],ylim,'k--');
            end
            axis tight
            if optDer
                title('Derivative and its radial quadrature','Interpreter','LaTex');
            else
                title('Function and its radial quadrature','Interpreter','LaTex');
            end
            
            xlabel('r');
        case 2
            t2 = linspace(0,1,fix(min(max(rho)*5,20000)));
            t1 = t2(and(t2>a/3,t2<b*crop));
            if optDer
                [~,err] = in.evalDer(t1);
            else
                [~,err] = in.eval(t1);
            end
            colors = get(gca,'ColorOrder');
            index  = get(gca,'ColorOrderIndex');
            loglog(t1,abs(err),'color',colors(index,:));
            hold on
            loglog(t1,t1*0+tol,'k--','HandleVisibility','off');
            loglog(ones(36,1)*a,10.^(linspace(-18,18,36)),'color',colors(index,:),'LineStyle',...
                '--','HandleVisibility','off');
            grid on
            axis tight
            ylim([1e-14 1]);
            title('Quadrature error','Interpreter','LaTex');
            xlabel('r','Interpreter','LaTex');
            set(gca,'ColorOrderIndex',mod(index,7)+1)
            
        case 3
            scatter(rho(2:end),abs(alpha0(2:end)).^2);
            hold on
            
            scatter(rho(2:end),abs(scal01).^2);
            
            set(gca,'yscale','log');
            
            grid on;
            axis tight
            if startFreq >0
                hold on
                plot([startFreq startFreq],ylim,'k--','LineWidth',2,'HandleVisibility','off');
            end
            title('Spectral power','Interpreter','LaTex');
            if isempty(scal01(isnan(scal01)))
                legend({'Optimal approximation','Truncature of Bessel-Fourier expansion'})
            end
            xlabel('$\rho$','Interpreter','LaTex');
        otherwise
            error('invalid argument n (value 1, 2, 3, or ''all'' accepted)')
    end
end

end


