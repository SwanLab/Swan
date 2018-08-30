function [xval,f0val,fval,kktnorm,xold2,outit] = mma(fun,x0,xmin,xmax,xold1,c,d,a0,a,options,varargin)
%---------------------------------------------------------------------
%  This is the file beammain.m.  Version September 2007.
%  Written by Krister Svanberg <krille@math.kth.se>.
%
%  This file contains a main program for using MMA to solve
%  a problem defined by the users files beaminit.m
%  (which must be run before beammain.m) and beam2.m.
%
%%%% If outeriter=0, the user should now calculate function values
%%%% and gradients of the objective- and constraint functions at xval.
%%%% The results should be put in f0val, df0dx, fval and dfdx:

outeriter = 0;
maxoutit = 1e4;
kkttol = 1e-6;
xval = x0;
xold2 = xold1;
if nargin >= 10
    if isfield(options,'outeriter')
        outeriter = options.outeriter;
    end
    if isfield(options,'maxoutit')
        maxoutit = options.maxoutit;
    end
    if isfield(options,'kkttol')
        kkttol = options.kkttol;
    end
end

outit = 0;
if outeriter < 0.5
    [f0val,df0dx,fval,dfdx] = fun(xval);
    low = xmin;
    upp = xmax;
    
     % DISPLAY INFO TO USER
%     display(sprintf('outit=%d, f0val=%.10f, fval=%.10f',outit,f0val,fval))
    if ~isempty(options.OutputFcn) % call output functions (update post_info, plots and save to file)
        optdata.compliance = f0val;
        optdata.gradient = df0dx;
        optdata.iter_ls = outit;
        optdata.incre_gamma = 0;
        optdata.lambda = zeros(length(c),1);
        for i = 1:length(options.OutputFcn)
            feval(options.OutputFcn{i},xval,optdata,varargin{:});
        end
    end
end
%
%%%% The iterations start:
m = length(c);
n = length(xval);
kktnorm = kkttol+10;
while kktnorm > kkttol && outit < maxoutit
    outit   = outit+1;
    outeriter = outeriter+1;
    %%%% The MMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
    %%%% Some vectors are updated:
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    %%%% The user should now calculate function values and gradients
    %%%% of the objective- and constraint functions at xval.
    %%%% The results should be put in f0val, df0dx, fval and dfdx.
    [f0val,df0dx,fval,dfdx] = fun(xval);
    %%%% The residual vector of the KKT conditions is calculated:
    [~,kktnorm] = ...
        kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    
    
    % DISPLAY INFO TO USER
%     display(sprintf('outit=%d, f0val=%.10f, fval=%.10f, kktnorm=%.10f',outit,f0val,fval,kktnorm))
    if ~isempty(options.OutputFcn) % call output functions (update post_info, plots and save to file)
        optdata.compliance = f0val;
        optdata.gradient = df0dx;
        optdata.iter_ls = outit;
        optdata.incre_gamma = kktnorm;
        optdata.lambda = lam;
        for i = 1:length(options.OutputFcn)
            feval(options.OutputFcn{i},xval,optdata,varargin{:});
        end
    end

end
%---------------------------------------------------------------------
