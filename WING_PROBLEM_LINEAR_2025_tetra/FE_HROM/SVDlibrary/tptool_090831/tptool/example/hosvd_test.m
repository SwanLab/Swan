clear;

% LPV model
lpv = {...
    @(p)1      @(p)p(1);
    @(p)p(2)^2 @(p)p(1)^2+p(2);
};
% number of states (size of the A matrix)
n = 4;
% parameter dependencies:
% dep(i,j,k) is 1 if Sp{i,j} depends on p(k)
dep = zeros([size(lpv) 2]);
dep(1,2,:) = [1 0];
dep(2,1,:) = [0 1];
dep(2,2,:) = [1 1];

% sampling intervals for each parameter
domain = [-5 5; -5 5];
% grid size: number of grid points for each parameter
gridsize = [11 11];

% hosvd_lpv
D1 = sampling_lpv(lpv, dep, domain, gridsize);
[S1 U1] = hosvd_lpv(D1, dep, gridsize);

% hosvd of a full array
f = @(p) [lpv{1,1}(p) lpv{1,2}(p); lpv{2,1}(p) lpv{2,2}(p)];
D2 = sampling_vec(f, domain, gridsize);
[S2 U2] = hosvd(D2, [1 1 0 0], 0.01);

% signs in U and S might be different
Sdiff = abs(S1(:)) - abs(S2(:));

% should be small
disp ' '
disp 'error (should be small):'
disp(norm(Sdiff))
