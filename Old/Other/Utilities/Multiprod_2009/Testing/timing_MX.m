function timing_MX
% TIMING_MX  Speed of MX as performed by MULTIPROD and by a nested loop.
%    TIMING_MX compares the speed of matrix expansion as performed by
%    MULTIPROD and an equivalent nested loop. The results are shown in the
%    manual (fig. 2).
%    Notice that MULTIPROD enables array expansion which generalizes matrix
%    expansion to arrays of any size, while the loop tested in this
%    function works only for this specific case, and would be much slower
%    if it were generalized to N-D arrays.


% Checking whether needed software exists
message = sysrequirements_for_testing('timeit');
if message
    disp ' ', error('testing_memory_usage:Missing_subfuncs', message)
end

% Matrix expansion example (fig. 2)
disp ' '
disp 'Timing matrix expansion (see MULTIPROD manual, figure 2)'
disp ' '

a = rand(2, 5);
b = rand(5, 3, 1000, 10);

fprintf ('Size of A:  %0.0fx%0.0f\n', size(a))
fprintf ('Size of B: (%0.0fx%0.0f)x%0.0fx%0.0f\n', size(b))
disp ' ', disp 'Please wait...'
disp ' '

f1 = @() loop(a,b);
f2 = @() multiprod(a,b);

t1 = timeit(f1)*1000;
fprintf('LOOP(A, B):      %10.4f milliseconds\n', t1)
t2 = timeit(f2)*1000;
fprintf('MULTIPROD(A, B): %10.4f milliseconds\n', t2)

disp ' '
fprintf('MULTIPROD performed matrix expansion %6.0f times faster than a plain loop\n', t1/t2)
disp ' '

function C = loop(A,B)
for i = 1:1000
    for j = 1:10
        C(:,:,i,j) = A * B(:,:,i,j);
    end
end