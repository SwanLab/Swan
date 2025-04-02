function testing_memory_usage
% TESTING_MEMORY_USAGE  Testing memory allocation by commands and engines.
%     This script requires MATLAB R2007a or higher (because it uses
%     BSXFUN), and at least 360 MB of free RAM.
%     Please open Task Manager (CTRL+ALT+DEL - "Processes" tag), to check
%     the memory used by the MATLAB process. Add a column showing "Peak
%     memory usage" (use menu "View" > "Select columns"). 
%     Before running this script, please close MATLAB to delete previous
%     value in the "Peak memory usage" column, then open it again.
%

% Paolo de Leva
% University of Rome, Foro Italico, Rome, Italy
% 2008 Dec 29

clear all; clc

% Checking whether needed software exists
message = sysrequirements_for_testing('bsxfun', 'genop', 'bsxtimes');
if message
    disp ' ', error('testing_memory_usage:Missing_subfuncs', message)
end

disp ' ', disp 'Read memory usage, then press any key'; pause

a1D = rand(1,1e6); 
disp ' '
fprintf ('Creating A1D. Size of A1D: %0.0f x %0.0f\n', size(a1D));
disp 'Read memory usage then press any key.'; pause

a2D = rand(100, 1e4); 
disp ' ' 
fprintf ('Creating A2D. Size of A2D: %0.0f x %0.0f\n', size(a2D));
disp 'Read memory usage then press any key.'; pause

% 1-D TRANSPOSE and RESHAPE require very little additional memory. T1D and
% A3D are just new pointers to the same data stack. See MATLAB help about
% memory usage. However, 2-D TRANSPOSE requires a new stack.
    t1D = a1D';
    disp ' '
    fprintf ('TRANSPOSE(A1D). Size of T1D: %0.0f x %0.0f\n', size(t1D));
    disp 'Read memory usage then press any key.'; pause

    t2D = a2D'; % Requires additional memory. A new stack.
    disp ' '
    fprintf ('TRANSPOSE(A2D). Size of T2D: %0.0f x %0.0f\n', size(t2D));
    disp 'Read memory usage then press any key.'; pause

    a3D = reshape(a2D, [1e4 1 100]);
    disp ' '
    fprintf ('RESHAPE(A2D, ...). Size of A3D: %0.0f x %0.0f x %0.0f\n', size(a3D));
    disp 'Read memory usage then press any key.'; pause

% A_PERM requires additional memory. A new stack.
    a_perm = permute(a2D, [2 1]); 
    disp ' '
    fprintf ('PERMUTE(A2D, ...). Size of A_PERM: %0.0f x %0.0f\n', size(a_perm));
    disp 'Read memory usage, then press any key.'; pause

 % A_BIG is a new stack, 5 times larger than A.
    a_big = a3D(:, ones(1,5), :); 
    disp ' '
    fprintf ('TONY''S TRICK. Size of A_BIG: %0.0f x %0.0f x %0.0f\n', size(a_big));
    disp 'Read memory usage, then press any key.'; pause
    
% B1 is a new stack as large as A_BIG. This statement does exactly the same as 
% B1 = (A3D(:,ones(1,5),:) .* A_BIG) but it does not require additional space 
% for expanding A3D, not even temporarily. However, BSXFUN is available
% only in MATLAB 2007a and later versions.
    b1 = bsxfun(@times, a3D, a_big); 
    disp ' '
    fprintf ('BSXFUN. Size of B1: %0.0f x %0.0f x %0.0f\n', size(b1));
    disp 'Read memory usage, then press any key.'; pause

% B2 is a new stack as large as A_BIG. This statement does not require
% additional space, and it is compatible with MATLAB 7 and earlier.
    b2 = bsx_times(a3D, a_big); % BSXFUN replacement by Doug Schwarz 
    disp ' '
    fprintf ('BSX_TIMES. Size of B2: %0.0f x %0.0f x %0.0f\n', size(b2));
    disp 'Read memory usage, then press any key.'; pause

% B3 is a new stack as large as A_BIG. This statement does exactly the same
% as the previous one, but it requires, temporarily, some additional space
% (perhaps for executing EVAL).
    b3 = genop(@times, a3D, a_big); % BSXFUN replacement by Doug Schwarz
    disp ' '
    fprintf ('GENOP. Size of B3: %0.0f x %0.0f x %0.0f\n', size(b3));
    disp 'Read memory usage, then press any key.'; pause

% B4 is a new stack as large as A_BIG, but the statement requires, temporarily,
% twice as much space in memory, to allocate both B3 and A3D(:,ones(1,5),:)
    b4 = (a3D(:,ones(1,5),:) .* a_big); 
    disp ' '
    fprintf ('T. TRICK EMBEDDED IN TIMES (.*). Size of B4: %0.0f x %0.0f x %0.0f\n', size(b4));
    disp 'Read memory usage, then press any key.'; pause

disp ' ', disp 'Thank you.'