function[MV,loc,K] = MVproduct(kernel_string,X,Y,tol,k,varargin)
% Inputs : 
% - kernel_string : one of the options listed in the switch loop
% - X, Y :          clouds of points
% - tol
% - k :             if a hankel kernel is selected in input, the wavenumber. Ignored if kernel_string is not
%                   a frequency dependent kernel.
% - varargin :      value pair 'lambda'. Enter values between 5 and 10 to
%                   rebalance the proportions of far and close interactions
% Outputs :     
% MV :              A handle such that f = MV(q) returns a fast approximation
%                   of f_i = \sum_{j} G(X_i - Y_j) q_j
% loc :             Sparse matrix containing the close interactions
%                   computed exactly
% K :               The kernel object corresponding to kernel_string


p = inputParser;
p.addOptional('lambda',1);
p.parse(varargin{:});
lambda = p.Results.lambda;
Nx = size(X,1);
Ny = size(Y,1);
a = lambda*1/sqrt(sqrt(Nx*Ny)); % Dichotomy instead ? 
switch(kernel_string)
    case '[log(r)]'
        K = LogKernel(1);
        [MV,~,loc] = offlineEBD(K,X,Y,a,tol);
    case '[H0(kr)]'
        K = H0Kernel(k);
        [MV,~,loc] = offlineEBD(K,X,Y,a,tol);
    case 'grady[log(r)]'
        K = LogKernel(1);
        [MVx,MVy,~,loc] = offline_dEBD(K,X,Y,a,tol);
        MV = {MVx,MVy};
    case 'grady[H0(kr)]'
        K = H0Kernel(k);
        [MVx,MVy,~,loc] = offline_dEBD(K,X,Y,a,tol);
        MV = {MVx,MVy};
end




end