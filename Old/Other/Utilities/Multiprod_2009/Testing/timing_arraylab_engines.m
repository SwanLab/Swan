function timing_arraylab_engines
% TIMING_ARRAYLAB_ENGINES  Testing for speed different ARRAYLAB engines.
% 
%    Testing different methods to multi-multiply A by B, which can be used
%    as kernels of the MULTIPROD function. The tested engines are specific 
%    for three input formats:
%
%    1)  A  is  N×(P×Q) (no block expansion is needed)
%        B  is  N×(Q×R)
%
%    2)  A  is  (P×Q)  (this array requires block expansion)
%        B  is  N×(Q×R)
%
%    3)  A  is  N×(P×Q)
%        B  is  (Q×R)  (this array requires block expansion)
%
%    The results obtained using two different personal computers are
%    reported in the manual of MULTIPROD 2.0 (Appendix A).

% Paolo de Leva  University of Rome, Foro Italico, Rome, Italy
% Jinhui Bai     Georgetown University, Washington, D.C.
% 2009 Gen 22

clear all

% Checking whether needed software exists
message = sysrequirements_for_testing('bsxfun','genop','bsxtimes','timeit');
if message
    disp ' ', error('timing_arraylab_engines:Missing_subfuncs', message)
end

disp ' ', 
disp '---------------------------- Experiment 1 ----------------------------'
N = 10000; P = 3; Q = 3; R = 1;  
[t(1,:), maxerror(1,:)] = timing(N,P,Q,R);

disp ' ', 
disp '---------------------------- Experiment 2 ----------------------------'
N = 1000; P = 3;  Q = 30; R = 1;  
[t(2,:), maxerror(2,:)] = timing(N,P,Q,R);

disp ' ', 
disp '---------------------------- Experiment 3 ----------------------------'
N = 1000; P = 9; Q = 10;  R = 3; 
[t(3,:), maxerror(3,:)] = timing(N,P,Q,R);

disp ' ', 
disp '---------------------------- Experiment 4 ----------------------------'
N = 100; P = 9; Q = 100;  R = 3; 
[t(4,:), maxerror(4,:)] = timing(N,P,Q,R);

disp ' ', 
disp '----------------------------   Summary   -----------------------------'
disp ' ', disp 'Execution time (ms)'
disp '     EXP 1     EXP 2     EXP 3     EXP 4'
disp (t'*1000)
disp 'Engine output error'
disp '     EXP 1     EXP 2     EXP 3     EXP 4'
fprintf('%10.2g%10.2g%10.2g%10.2g\n', maxerror')
disp 'NOTE: Engine output error should be zero or very close to EPS.'
fprintf ('      EPS ='); disp (eps)



function [t, maxerror] = timing(N,P,Q,R)

a0 = rand(1, P, Q); 
b0 = rand(1, Q, R);
a = a0(ones(1,N),:,:); % Cloning along first dimension
b = b0(ones(1,N),:,:); 
[n1 p q1] = size(a); % Reads third dim even if it is 1
[n2 q2 r] = size(b); %
disp ' '
disp        'Array  Size      Size               Number of elements'
fprintf (1, 'A      Nx(PxQ)   %0.0f x (%0.0f x %0.0f) %8.0f\n', [n1 p q1  numel(a)])
fprintf (1, 'B      Nx(QxR)   %0.0f x (%0.0f x %0.0f) %8.0f\n', [n2 q2 r  numel(b)])
disp ' '

% Standard loops (Probably they exploit MATLAB JIT accelerator)
disp 'ARRAYLAB by means of plain LOOPS'
f1  = @() loop1(a,b); 
f1b = @() loop1b(a,b);
f2  = @() loop2(a,b);
t = [timeit(f1), timeit(f1b), timeit(f2)];
disp (t)

% Method used in MULTIPROD 1.3, which uses "cloning", .* and SUM
% ("Cloning" means "singleton expansion"; see also BSXFUN help)
disp 'MULTIPROD VERSION 1 (column-by-column cloning) AND 2 (BSXFUN engine)'
disp '       1.3       1.3      1.31      1.32      1.33       2.1'
f3 = @() multiprod13(a, b, [2 3]);
f4 = @()  arraylab13(a, b, 2, 3); 
f5 = @() arraylab131(a, b, 2, 3); 
f6 = @() arraylab132(a, b, 2, 3); 
f7 = @() arraylab133(a, b, 2, 3); 
f8 = @()   multiprod(a, b, [2 3]);

lenT = length(t);
t = [t, timeit(f3),timeit(f4),timeit(f5),timeit(f6),timeit(f7),timeit(f8)];
disp (t(lenT+1:end))

% (PERMUTE) --> 4D RESHAPE --> SX --> TIMES --> SUM
disp 'ARRAYLAB by means of (PERMUTE) --> 4D-RESHAPE --> SX --> TIMES --> SUM'
disp ' TONY''S T.     GENOP  BSXTIMES    BSXFUN   BSXFUN2   BSXFUN3'
f9  = @() resh4D_clone(a, b);    % SX with Tony's trick
f10 = @() resh4D_genop(a, b);    % SX+TIMES with GENOP  (BSXFUN replacement)
f11 = @() resh4D_bsxtimes(a, b); % SX+TIMES with BSXTIMES (BSXFUN replacement)
f12 = @() resh4D_bsxfun(a, b);   % SX+TIMES with BSXFUN
f13 = @() permA_resh4D_bsxfun(a, b);
f14 = @() permAB_resh4D_bsxfun(a, b);
lenT = length(t);
t = [t, timeit(f9),timeit(f10),timeit(f11),timeit(f12),timeit(f13),timeit(f14)];
disp (t(lenT+1:end))

% Techniques for virtual matrix expansion (MX)
disp 'MATRIX EXPANSION (MX) TECHNIQUES (MULTIPROD 2 USES 2D-RESHAPE)'
disp '                          LOOP   4D-RESHAPE   2D-RESHAPE  MULTIPROD 2'
a2 = shiftdim(a(1,:,:)); % Single matrix
b2 = shiftdim(b(1,:,:)); % Single matrix
f15 = @()         loop2_expA(a2, b);     % PLAIN LOOP2
f16 = @() resh4D_bsxfun_expA(a2, b);     % 4D RESH --> SX --> TIMES --> SUM
f17 = @()     squashB_mtimes(a2, b);     % PERMUTE --> 2D RESH --> MTIMES
f18 = @() multiprod(a2, b, [1 2],[2 3]); % Generalized SQUASH_MTIMES
f19 = @()         loop2_expB(a, b2);
f20 = @() resh4D_bsxfun_expB(a, b2);
f21 = @()     squashA_mtimes(a, b2);
f22 = @() multiprod(a, b2, [2 3],[1 2]);
lenT = length(t);
t = [t, timeit(f15),timeit(f16),timeit(f17),timeit(f18)];
fprintf('MATRIX * 3D ARRAY')
fprintf(1, '%13.4f', t(lenT+1:end))
disp ' '
lenT = length(t);
t = [t, timeit(f19),timeit(f20),timeit(f21),timeit(f22)];
fprintf('3D ARRAY * MATRIX')
fprintf(1, '%13.4f', t(lenT+1:end))
disp ' '

% Comparing function output (differences)
c0 = f1();
c = c0 - f1b(); maxerror(1)  = max(c(:));
c = c0 - f2();  maxerror(2)  = max(c(:));
c = c0 - f3();  maxerror(3)  = max(c(:));
c = c0 - f4();  maxerror(4)  = max(c(:));
c = c0 - f5();  maxerror(5)  = max(c(:));
c = c0 - f6();  maxerror(6)  = max(c(:));
c = c0 - f7();  maxerror(7)  = max(c(:));
c = c0 - f8();  maxerror(8)  = max(c(:));
c = c0 - f9();  maxerror(9)  = max(c(:));
c = c0 - f10(); maxerror(10) = max(c(:));
c = c0 - f11(); maxerror(11) = max(c(:));
c = c0 - f12(); maxerror(12) = max(c(:));
c = c0 - f13(); maxerror(13) = max(c(:));
c = c0 - f14(); maxerror(14) = max(c(:));
c = c0 - f15(); maxerror(15) = max(c(:));
c = c0 - f16(); maxerror(16) = max(c(:));
c = c0 - f17(); maxerror(17) = max(c(:));
c = c0 - f18(); maxerror(18) = max(c(:));
c = c0 - f19(); maxerror(19) = max(c(:));
c = c0 - f20(); maxerror(20) = max(c(:));
c = c0 - f21(); maxerror(21) = max(c(:));
c = c0 - f22(); maxerror(22) = max(c(:));


function c = loop1(a, b)
% WARNING: You can use LOOP2 without PERMUTE if  A  is  (P×Q)×N  and
%                                                B  is  (Q×R)×N
[n p q] = size(a);
r = size(b, 3);
c = zeros([n p r]);
for i = 1 : n
    c(i,:,:)= reshape(a(i,:,:), p, q)...
              * ...
              reshape(b(i,:,:), q, r);
end

function c = loop1b(a, b)
% Same as LOOP 1 but commands are not nested
% (I have read somewhere that JIT does not work with nested commands)
[n p q] = size(a);
r = size(b, 3);
c = zeros([n p r]);
for i = 1 : n
    a2 = reshape(a(i,:,:), p, q);
    b2 = reshape(b(i,:,:), q, r);
    c(i,:,:)= a2 * b2;
end

function c = loop2(a, b)
% Warning:  PERMUTE is not needed if  A  is  (P×Q)×N  and
%                                     B  is  (Q×R)×N
[n p q] = size(a);
r = size(b, 3);
a = permute(a, [2 3 1]);
b = permute(b, [2 3 1]);
    c = zeros([p, r, n]);
    for i = 1 : n
        c(:,:,i) = a(:,:,i) * b(:,:,i);
    end
c = permute(c, [3 1 2]);
     
function c = resh4D_clone(a, b)
% By Paolo de Leva 
% Based on RESH4D_BSXFUN. Uses "Tony's trick" instead of BSXFUN. Tony's
% trick is an index vectorization trick used in REPMAT and MULTIPROD 1.3 to
% perform singleton expansion)
% WARNING: Tony's trick is memory expensive.
    [n p q] = size(a);
    [n q r] = size(b);
    b = reshape(b, [n 1 q r]);
    indexP = ones(1, p); % "Cloned" indexes
    indexR = ones(1, r);
    c = sum(a(:,:,:,indexR) .* b(:,indexP,:,:), 3);
    c = reshape(c, [n p r]);
    
function c = resh4D_bsxfun(a, b)
% By Jinhui Bai
    [n p q] = size(a);
    [n q r] = size(b);
    b = reshape(b, [n 1 q r]);
    c = sum(bsxfun(@times,a,b), 3);
    c = reshape(c, [n p r]);
    
function c = permA_resh4D_bsxfun(a, b)
% This version of RESH4D_BSXFUN permutes only A to make it possible a SUM
% along 2nd dimension (supposedly faster than along 3rd dimension). 
% By Jinhui Bai
    [n p q] = size(a);
    [n q r] = size(b);
    a = permute(a, [1 3 2]); % N x Q x P
    b = reshape(b, [n q 1 r]);
    c = sum(bsxfun(@times,a,b), 2);
    c = reshape(c, [n p r]);

function c = permAB_resh4D_bsxfun(a, b)
% This version of RESH4D_BSXFUN permutes both A and B to make it possible
% a SUM along 1st dimension (supposedly faster than along higher dim.). 
% By Jinhui Bai
    [n p q] = size(a);
    [n q r] = size(b);
    a = permute(a, [3 2 1]); % Q x P x N
    b = permute(b, [2 3 1]); % Q x R x N
        a = reshape(a, [q p 1 n]);
        b = reshape(b, [q 1 r n]);
        c = sum(bsxfun(@times,a,b), 1);
        c = reshape(c, [p r n]);
    c = permute(c, [3 1 2]); % N x P x R
    
function c = resh4D_genop(a, b)
% Based on RESH4D_BSXFUN. Uses GENOP by Doug Schwarz (MATLAB Central file
% #10333) as a replacement of BSXFUN.
    [n p q] = size(a);
    [n q r] = size(b);
    b = reshape(b, [n 1 q r]);
    c = sum(genop(@times,a,b), 3);
    c = reshape(c, [n p r]);
    
function c = resh4D_bsxtimes(a, b)
% Based on RESH4D_BSXFUN. Uses a BSXFUN substitute by Douglas Schwarz
% (MATLAB Central file #23005).
    [n p q] = size(a);
    [n q r] = size(b);
    b = reshape(b, [n 1 q r]);
    c = sum(bsx_times(a,b), 3);
    c = reshape(c, [n p r]);

function c = loop2_expA(a, b)
% Virtual single-block expansion
      p = size(a, 1);
[n q r] = size(b);
b = permute(b, [2 3 1]);
    c = zeros([p, r, n]);
    for i = 1 : n
        c(:,:,i) = a * b(:,:,i);
    end
c = permute(c, [3 1 2]);

function c = loop2_expB(a, b)
% Virtual single-block expansion
[n p q] = size(a);
      r = size(b, 2);
a = permute(a, [2 3 1]);
    c = zeros([p, r, n]);
    for i = 1 : n
        c(:,:,i) = a(:,:,i) * b;
    end
c = permute(c, [3 1 2]);

function c = resh4D_bsxfun_expA(a, b)
% Virtual single-block expansion
% Based on RESH4D_BSXFUN (singleton expansion).
% by Paolo de Leva
      [p q] = size(a);
    [n q r] = size(b);
    a = reshape(a, [1 p q 1]);
    b = reshape(b, [n 1 q r]);
    c = sum(bsxfun(@times,a,b), 3);
    c = reshape(c, [n p r]);
    
function c = resh4D_bsxfun_expB(a, b)
% Virtual single-block expansion 
% Based on RESH4D_BSXFUN (singleton expansion).
% by Paolo de Leva
    [n p q] = size(a);
      [q r] = size(b);
    b = reshape(b, [1 1 q r]);
    c = sum(bsxfun(@times,a,b), 3);
    c = reshape(c, [n p r]);   

function c = squashB_mtimes(a, b)
% Virtual single-block expansion
% Instead of expanding both A and B, it just squashes B from 3-D to 2-D. It
% avoids using additional memory for expansion. Based on code suggested by
% Jinhui Bai to deal with this format:
%     A  is  P×Q
%     B  is (Q×R)×N
%
% NOTE : B = PERMUTE(B, [2 1 3]), from N×(Q×R) to Q×N×R has been shown to 
%        be faster than B = PERMUTE(B, [2 3 1]), from N×(Q×R) to (Q×R)×N.
%        See TIMING_MATLAB_COMMANDS.M.
% By Paolo de Leva
b = permute(b, [2 1 3]); % Q×N×R
    p = size(a,1);
    [q n r] = size(b);
    b = reshape(b, [q, n*r]); % squashing B (equivalent to B(:,:))
    c = reshape(a*b, [p n r]); 
c = permute(c, [2 1 3]); % N×(P×R)
    
function c = squashA_mtimes(a, b)
% Virtual single-block expansion
% Instead of expanding both A and B, it just squashes A from 3-D to 2-D. It
% avoids using additional memory for expansion. Based on code suggested by
% Jinhui Bai to deal with this format:
%     A  is (P×Q)×N
%     B  is  Q×R
% By Paolo de Leva
    [n p q] = size(a);
    r = size(b, 2);
    a = reshape(a, [n*p, q]); % squashing A
    c = reshape(a*b, [n p r]);
