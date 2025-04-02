classdef Quad2D
    % Object that represents the approximation of a function in the form
    % f(x) \approx offset + \sum_{\nu = 1}^{N_{\xi}} \hat{w}_{\nu}e^{i x \cdot \xi_{\nu}}
    % valid for a < |x| < 1, where \xi_{\nu} is a list of quadrature points
    % and \hat{w}_{\nu} a list of quadrature weights.
    properties (Constant,Hidden)
        % Constants developped to help choosing the number of points for
        % the circular discretisation of the spectrum.
        gamma = exp(1)/2*0.8; % In the article, I show that exp(1)/2 is sufficient. 
        % In practice, exp(1)/2*0.8 seems to be enough. 
    end
    properties (Access = public)
        %rq; % The radial quadrature
        w_nu; % Quadrature weights
        xi_nu; % Quadrature points
        time; % The time it took to compute this object
        offset; % The offset as explained in the description 
        % (weight associated to 0 frequency). 
        tol;
    end
    
    methods
        function[q2d] = Quad2D(radialQuad)
            if nargin == 0
                q2d.w_nu = [];
                q2d.xi_nu = [];
                q2d.time = 0;
                q2d.offset = 0;
                q2d.tol = 0;
            else
                % Constructor
                ttime = tic;
                alpha= radialQuad.alpha;
                rho = radialQuad.rho;
                P = length(rho);
                assert(rho(1)==0);
                tol = radialQuad.tol;
                q2d.offset = alpha(1);
                alpha = alpha(2:end);
                rho = rho(2:end);
                Ns = fix((Quad2D.gamma * rho  + (4*log(P*abs(alpha)/(tol+10^(-10)))))/2 + 1)*2+1;
                N = sum(Ns);
                xxi_nu = zeros(N,2);
                ww_nu = zeros(N,1);
                i1 = 1;
                % Loop on frequencies and add quadrature points as discretized
                % circle of radius this frequency.
                for p = 1:length(rho)
                    i2 = Ns(p) + i1 - 1;
                    theta = (1:Ns(p))' * (2*pi)/Ns(p) ;
                    xxi_nu(i1:i2,1) = rho(p)*cos(theta);
                    xxi_nu(i1:i2,2) = rho(p)*sin(theta);
                    ww_nu(i1:i2) = alpha(p)/Ns(p) * ones(Ns(p),1);
                    i1 = i2 + 1;
                end
                q2d.w_nu = ww_nu;
                q2d.xi_nu = xxi_nu;
                q2d.time = toc(ttime);
                q2d.tol = tol;
            end
            
        end
        function[] = disp(this)
            fprintf('2D quadrature \n')
            fprintf('Number of quadrature points : %s\n',num2str(length(this.w_nu)))
        end
        % Number of frequency samples
        function[out] = Nxi(this)
            out = length(this.w_nu);
        end
        % Value of the functions at X(i,:) - Y(j,:)
        function[el] = elem(this,X,Y,i,j)
            qj = zeros(size(X,1),1);
            qj(j) = 1;
            col = this.conv(X,Y,qj);
            el = col(i);
        end
        function[] = show(this)
            figure
            subplot(1,3,1);
            scatter3(this.xi_nu(:,1),this.xi_nu(:,2),log(abs(this.w_nu)));
            xlabel('\xi_x');  
            ylabel('\xi_y');
            zlabel('log(|\omega_\nu|)')
            title('Quadrature coefficients (absolute values)')
            subplot(1,3,2);
            x = linspace(-1,1,100);
            y = linspace(-1,1,100);
            [X, Y] = meshgrid(x,y);
            Z = this.value([X(:) Y(:)]);
            Z = reshape(Z,size(X,1),size(X,2));
            surf(X,Y,real(Z));
            title('real part');
            subplot(1,3,3);
            surf(X,Y,imag(Z));
            title('Imaginary part');
        end
        % Addition of two 2D quadratures
        function[c] = plus(a,b)
            t = tic;
            c = Quad2D;
            xi_complex_a = a.xi_nu(:,1) + 1i*a.xi_nu(:,2);
            xi_complex_b = b.xi_nu(:,1) + 1i*b.xi_nu(:,2);
            [xi_complex_c,c.w_nu] = mergeQuads(xi_complex_a,a.w_nu,xi_complex_b,b.w_nu);
            c.xi_nu = [real(xi_complex_c) imag(xi_complex_c)];            
            c.offset = a.offset + b.offset;
            c.tol = a.tol + b.tol;
            c.time = toc(t);
        end
        % Value at a point in R^2
        function[out] = value(this,x)
            xxi_nu = this.xi_nu;
            ww_nu = this.w_nu;
            Nxi = length(xxi_nu);
            ttol = this.tol/100 + 1e-10;
            if Nxi >0
                out = nufft2d3(Nxi, xxi_nu(:,1), xxi_nu(:,2), ...
                    ww_nu, +1, ttol, size(x,1), x(:,1), x(:,2)) + this.offset;
            else
                out = 0*size(x,1);
            end
        end
        function[dx,dy] = grad(this)
            dx = this;
            dx.w_nu = this.w_nu.*(1i*this.xi_nu(:,1));
            dx.offset = 0;
            dy = this;
            dy.w_nu = this.w_nu.*(1i*this.xi_nu(:,2));
            dy.offset = 0;
        end
        % Fast convolution
        function[q,time] = conv(this,x,y,V)
            % This function computes \sum_{j=1}^N f(x(i)-y(j)) V(j)
            % where f is the function approximated by this object.
            % This computes only accurately the far contribution, i.e. |x(i) - y(j)| > a
            % while the close interactions have to be corrected.
            % |x(i) - y(j)| has to be <= 1
            time = tic;
            xxi_nu = this.xi_nu;
            ww_nu = this.w_nu;
            Nxi = length(xxi_nu);
            ttol = min(this.tol/10 + 1e-10,1e-6);
            if Nxi > 0
                %% Space to Fourier
                V_nu = nufft2d3(size(y,1), y(:,1), y(:,2), ...
                    V, -1, ttol, Nxi, xxi_nu(:,1), xxi_nu(:,2) );
                %% Convolution becomes multiplication of Fourier weights
                fV_nu = V_nu.*ww_nu;
                %% Back from Fourier to Space
                q = nufft2d3(Nxi, xxi_nu(:,1), xxi_nu(:,2), ...
                    fV_nu, +1, ttol, size(x,1), x(:,1), x(:,2));
                time = toc(time);
                %% Adding offset
                q = q + this.offset*sum(V);
            else
                q = 0*size(x,1);
                time = toc(time);
            end
        end
        
    end
end

