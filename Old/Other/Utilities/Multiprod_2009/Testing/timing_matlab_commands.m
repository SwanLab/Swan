function timing_matlab_commands
% TIMING_MATLAB_COMMANDS  Testing for speed different MATLAB commands.
% 
% Main conclusion: RESHAPE and * (i.e. MTIMES) are very quick!

% Paolo de Leva
% University of Rome, Foro Italico, Rome, Italy
% 2008 Dec 24

clear all

% Checking whether needed software exists
if ~exist('bsxfun', 'builtin')
    message = sysrequirements_for_testing('bsxmex', 'timeit');
else
    message = sysrequirements_for_testing('timeit');
end
if message
    disp ' ', error('timing_matlab_commands:Missing_subfuncs', message)
end

disp ' '
disp '---------------------------------- Experiment 1 ----------------------------------'
N = 10000; P = 3; Q = 3; R = 1;  
timing(N,P,Q,R);

disp '---------------------------------- Experiment 2 ----------------------------------'
N = 1000; P = 3;  Q = 30; R = 1;  
timing(N,P,Q,R);

disp '---------------------------------- Experiment 3 ----------------------------------'
N = 1000; P = 9; Q = 10;  R = 3; 
timing(N,P,Q,R);

disp '---------------------------------- Experiment 4 ----------------------------------'
N = 100; P = 9; Q = 100;  R = 3; 
timing(N,P,Q,R);

disp '---------------------------------- Experiment 5 ----------------------------------'
disp ' '
timing2(4, 10000);
timing2(200, 200);
timing2(10000, 4);

disp '---------------------------- Experiment 6 ----------------------------'
disp ' '
a = rand(4096, 4096);
fprintf ('Size of A:  %0.0f x %0.0f\n', size(a))
disp ' '
disp '   SUM(A,1)  SUM(A,2)'
f1 = @() sum(a, 1);
f2 = @() sum(a, 2);
disp ([timeit(f1), timeit(f2)])

clear all
b = rand(256, 256, 256);
fprintf ('Size of B:  %0.0f x %0.0f x %0.0f\n', size(b))
disp ' '
disp '   SUM(B,1)  SUM(B,2)  SUM(B,3)'
f1 = @() sum(b, 1);
f2 = @() sum(b, 2);
f3 = @() sum(b, 3);
disp ([timeit(f1), timeit(f2), timeit(f3)])

disp '---------------------------- Experiment 7 ----------------------------'
disp ' '
a = rand(101,102,103);
fprintf ('Size of A:  %0.0f x %0.0f x %0.0f\n', size(a))
disp ' '
disp 'Moving last dimension to first dimension:'
disp 'PERMUTE(A,[3 2 1])  PERMUTE(A,[3 1 2])  SHIFTDIM(A,2)'
disp '(SWAPPING)          (SHIFTING)          (SHIFTING)'
f1 = @() permute(a, [3 2 1]);
f2 = @() permute(a, [3 1 2]);
f3 = @() shiftdim(a, 2);
fprintf(1, '%8.2g            ', [timeit(f1), timeit(f2), timeit(f3)])
disp ' ', disp ' '
a2 = f1(); s = size(a2);
a2 = f2(); s(2,:) = size(a2);
a2 = f3(); s(3,:) = size(a2);
disp (s)

disp 'Moving first dimension to last dimension:'
disp 'PERMUTE(A,[3 2 1])  PERMUTE(A,[2 3 1])  SHIFTDIM(A,1)'
disp '(SWAPPING)          (SHIFTING)          (SHIFTING)'
f1 = @() permute(a, [3 2 1]);
f2 = @() permute(a, [2 3 1]);
f3 = @() shiftdim(a, 1);
fprintf(1, '%8.2g            ', [timeit(f1), timeit(f2), timeit(f3)])
disp ' ', disp ' '
a2 = f1(); s = size(a2);
a2 = f2(); s(2,:) = size(a2);
a2 = f3(); s(3,:) = size(a2);
disp (s)

disp ' '
a = rand(21,22,23,24,25);
fprintf ('Size of A:  %0.0f x %0.0f x %0.0f x %0.0f x %0.0f\n', size(a))
disp ' '
disp 'Moving 4th dimension to 1st dimension:'
disp 'PERMUTE(A,[4 2 3 1 5])  PERMUTE(A,[4 1 2 3 5])  PERMUTE(A,[4 5 1 2 3])'
disp '(SWAPPING)              (PARTIAL SHIFTING)      (SHIFTING)'
f1 = @() permute(a, [4 2 3 1 5]);
f2 = @() permute(a, [4 1 2 3 5]);
f3 = @() permute(a, [4 5 1 2 3]);
fprintf(1, '%8.2g                ', [timeit(f1), timeit(f2), timeit(f3)])
disp ' ', disp ' '
a2 = f1(); s = size(a2);
a2 = f2(); s(2,:) = size(a2);
a2 = f3(); s(3,:) = size(a2);
disp (s)

disp 'Moving 2nd dimension to 5th dimension:'
disp 'PERMUTE(A,[1 5 3 4 2])  PERMUTE(A,[1 3 4 5 2])  PERMUTE(A,[3 4 5 1 2])'
disp '(SWAPPING)              (PARTIAL SHIFTING)      (SHIFTING)'
f1 = @() permute(a, [1 5 3 4 2]);
f2 = @() permute(a, [1 3 4 5 2]);
f3 = @() permute(a, [3 4 5 1 2]);
fprintf(1, '%8.2g                ', [timeit(f1), timeit(f2), timeit(f3)])
disp ' ', disp ' '
a2 = f1(); s = size(a2);
a2 = f2(); s(2,:) = size(a2);
a2 = f3(); s(3,:) = size(a2);
disp (s)

disp '---------------------------- Experiment 8 ----------------------------'
disp ' '
a =rand(101,102,103);
order = [1 2 3];
shape = [101,102,103];
f1 = @() perm(a,order);
f2 = @() ifpermute(a,order);
f3 = @() ifpermute2(a,order);
f4 = @() resh(a,shape);
f5 = @() ifreshape(a,shape);
f6 = @() ifreshape2(a,shape);
disp 'COMPARING STATEMENTS THAT DO NOTHING!'
disp ' '
fprintf ('Size of A:  %0.0f x %0.0f x %0.0f\n', size(a))
disp ' '
disp 'ORDER = [1 2 3]       % (keeping same order)'
disp 'SHAPE = [101,102,103] % (keeping same shape)'
disp ' '
fprintf (1,'PERMUTE(A,ORDER) ..........................................  %0.4g\n', timeit(f1))
fprintf (1,'IF ~ISEQUAL(ORDER,1:LENGTH(ORDER)), A=PERMUTE(A,ORDER); END  %0.4g\n', timeit(f2))
fprintf (1,'IF ~ISEQUAL(ORDER,1:3),             A=PERMUTE(A,ORDER); END  %0.4g\n', timeit(f3))
disp ' '
fprintf (1,'RESHAPE(A,SHAPE) ..........................................  %0.4g\n', timeit(f4))
fprintf (1,'IF ~ISEQUAL(SHAPE,SIZE(A)), A=RESHAPE(A,SHAPE); END .......  %0.4g\n', timeit(f5))
fprintf (1,'IF ~ISEQUAL(SHAPE,SHAPE),   A=RESHAPE(A,SHAPE); END .......  %0.4g\n', timeit(f5))
disp ' '


function a=perm(a, order)
a=permute(a, order);
function a=resh(a,shape)
a=reshape(a,shape);
function a=ifpermute(a, order)
if ~isequal(order, 1:length(order)), a=permute(a,order); end
function a=ifreshape(a, shape)
if ~isequal(shape, size(a)), a=reshape(a,shape); end
function a=ifpermute2(a, order)
if ~isequal(order, 1:3), a=permute(a,order); end
function a=ifreshape2(a, shape)
if ~isequal(shape, shape), a=reshape(a,shape); end


function timing(N,P,Q,R)

a0 = rand(1, P, Q); 
b0 = rand(1, Q, R);
a = a0(ones(1,N),:,:); % Cloning along first dimension
b = b0(ones(1,N),:,:); % Cloning along first dimension
[n1 p q1] = size(a); % reads third dim even if it is 1.
[n2 q2 r] = size(b); % reads third dim even if it is 1.
disp ' '
disp        'Array  Size      Size               Number of elements'
fprintf (1, 'A      Nx(PxQ)   %0.0f x (%0.0f x %0.0f) %8.0f\n', [n1 p q1  numel(a)])
fprintf (1, 'B      Nx(QxR)   %0.0f x (%0.0f x %0.0f) %8.0f\n', [n2 q2 r  numel(b)])
f1 = @() permute(a, [2 3 1]);
f2 = @() permute(a, [1 3 2]);
f3 = @() permute(a, [2 1 3]);
f4 = @() permute(a, [1 2 3]);
f5 = @() permute(b, [2 3 1]);
f6 = @() permute(b, [1 3 2]);
f7 = @() permute(b, [2 1 3]);
f8 = @() permute(b, [1 2 3]);
disp ' '
disp '   PERMUTE(A,[2 3 1])  PERMUTE(A,[1 3 2])  PERMUTE(A,[2 1 3])  PERMUTE(A,[1 2 3])'
fprintf(1, '%20.5f', [timeit(f1), timeit(f2), timeit(f3), timeit(f4)])
disp ' '
disp '   PERMUTE(B,[2 3 1])  PERMUTE(B,[1 3 2])  PERMUTE(B,[2 1 3])  PERMUTE(B,[1 2 3])'
fprintf(1, '%20.5f', [timeit(f5), timeit(f6), timeit(f7), timeit(f8)])
disp ' '
disp ' '
disp '   RESHAPE(A,[N*P Q])  RESHAPE(B,[N R Q])  RESHAPE(B,[N 1 R Q])'
f1 = @() reshape(a, [N*P Q]);
f2 = @() reshape(b, [N R Q]);
f3 = @() reshape(b, [N 1 R Q]);
fprintf(1, '%20.5f', [timeit(f1), timeit(f2), timeit(f3)])
disp ' '
f1 = @() a .* a;
f2 = @() bsxfun(@times, a, a);
f3 = @() b .* b;
f4 = @() bsxfun(@times, b, b);
disp ' '
disp '              A .* A   BSXFUN(@TIMES,A,A)'
fprintf(1, '%20.5f%20.5f\n', [timeit(f1), timeit(f2)])
disp '              B .* B   BSXFUN(@TIMES,B,B)'
fprintf(1, '%20.5f%20.5f\n', [timeit(f3), timeit(f4)])
if R==1
    disp ' '
    disp '   NOTE: If R=1 then RESHAPE(B,[N R Q]) is equivalent to'
    disp '                     PERMUTE(B,[1 3 2]) but much faster!'
    disp '                     (at least on my system)'
end
disp ' '


function timing2(P,Q)

a = rand(P, Q); 
b = rand(Q, 1);
fprintf ('Size of A:  %0.0f x %0.0f\n', size(a))
fprintf ('Size of B:  %0.0f x %0.0f\n', size(b))
disp ' '
disp '        A * B   TONY''S TRICK     BSXFUN'
f1 = @() a * b;
f2 = @() clone_multiply_sum(a, b', P);
f3 = @() sum(bsxfun(@times, a, b'), 2);
fprintf(1, '%13.5f', [timeit(f1), timeit(f2), timeit(f3)])
disp ' '
disp ' '
c = f1() - f2(); 
d = max(c(:)); 
if  d > eps*20
    disp 'There is an unexpected output difference:';
    disp (d);
end

function c = clone_multiply_sum(a,b,P)
c = sum(a .* b(ones(1,P),:), 2);