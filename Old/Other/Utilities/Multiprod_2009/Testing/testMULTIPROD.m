function testMULTIPROD
% TESTMULTIPROD  Testing function MULTIPROD
%    To perform a series of tests of function MULTIPROD, use
%    TESTMULTIPROD (without arguments)

global numberofchecks  min_of_all_e  max_of_all_e  syntaxes

format compact;
min_of_all_e = 0;
max_of_all_e = 0;
Nsyn = 11; % Number of MULTIPROD syntaxes
syntaxes = '     1    1s   1xs    2a    2b    3a    3b   4ab   4abs   4b    4c';
numberofchecks = zeros(1, Nsyn);

% Initial size of matrices
% (Q = R required, unless P=Q=1 or R=S=1)
p=6; q=3;
r=3; s=4;

% 0 singletons
    check ([p q], [r s]);

% 1 singleton
      check ([1 q], [r s]);
    % check ([p 1], [r s]); % Not allowed
    % check ([p q], [1 s]); % Not allowed
      check ([p q], [r 1]);

% 2 singletons
      check ([1 1], [r s]);
      check ([p q], [1 1]);
      check ([1 q], [r 1]); % inner product
      check ([p 1], [1 s]); % outer product
    % check ([1 q], [1 s]); % Not allowed
    % check ([p 1], [r 1]); % Not allowed

% 3 singletons
    check ([1 1], [1 s]);
    check ([1 1], [r 1]);
    check ([1 q], [1 1]);
    check ([p 1], [1 1]);

% 4 singletons
    check ([1 1], [1 1]);

% EMPTY MATRICES    
%-------------------------------------------

% 0 singletons
    check ([0 q], [r s]);
    check ([p q], [r 0]);
    check ([p 0], [0 s]);
    check ([0 0], [0 s]);
    check ([p 0], [0 0]);
    check ([0 0], [0 0]);

% 1 singleton
      check ([1 q], [r 0]);
      check ([1 0], [0 s]);
      check ([1 0], [0 0]);
    % check ([p 1], [r s]); % Not allowed
    % check ([p q], [1 s]); % Not allowed
      check ([0 q], [r 1]);
      check ([p 0], [0 1]);
      check ([0 0], [0 1]);

% 2 singletons
      check ([1 1], [0 s]);
      check ([1 1], [r 0]);
      check ([1 1], [0 0]);

      check ([0 q], [1 1]);
      check ([p 0], [1 1]);
      check ([0 0], [1 1]);

      check ([1 0], [0 1]); % inner product
      check ([0 1], [1 s]); % outer product
      check ([p 1], [1 0]); % outer product
      check ([0 1], [1 0]); % outer product
    % check ([1 q], [1 s]); % Not allowed
    % check ([p 1], [r 1]); % Not allowed

% 3 singletons
    check ([1 1], [1 0]);
    check ([1 1], [0 1]);
    check ([1 0], [1 1]);
    check ([0 1], [1 1]);

disp ' '
disp '------------------------------------------------------------------'
disp ' '
fprintf 'Total number of checked multiproducts:', disp (sum(numberofchecks))
disp ' ', disp 'Number of checked multiproducts for each syntax:'
disp (syntaxes)
str = ('%6.0f')';
str = str(:, ones(1,Nsyn));
formatstr = [str(:)' '\n'];
fprintf(1, formatstr, numberofchecks);
disp ' '
disp ( ['Minimun error:    ' num2str(min_of_all_e)] )
disp ( ['Maximum error:    ' num2str(max_of_all_e)] )
disp ( ['MATLAB precision: ' num2str(eps)] )


%--------------------------------------------------------------------------
function check(blsizeA, blsizeB)

if     all(blsizeA==1), blsizeC = blsizeB; % scalars in A
elseif all(blsizeB==1), blsizeC = blsizeA; % scalars in B
else                    blsizeC = [blsizeA(1) blsizeB(2)]; % No scalars
end

a = zeros([blsizeA 5 2]);
b = zeros([blsizeB 5 2]);
c = zeros([blsizeC 5 2]);
for d1 = 1:5
    for d2 = 1:2 
        a0 = rand( blsizeA ); a(:,:,d1,d2) = a0; 
        b0 = rand( blsizeB ); b(:,:,d1,d2) = b0; 
        c0 = a0 * b0;         c(:,:,d1,d2) = c0;
    end
end
a1  = a(:,:,:,:); b1  = b(:,:,:,:); c1  = c; % 4-D by 4-D
a2  = a(:,:,1,:); b2  = b(:,:,:,:); c2  = c; % 4-D by 4-D
a3  = a(:,:,:,1); b3  = b(:,:,:,:); c3  = c; % 3-D by 4-D
a4  = a(:,:,1,1); b4  = b(:,:,:,:); c4  = c; % 2-D by 4-D
a5  = a(:,:,:,:); b5  = b(:,:,1,:); c5  = c; % 4-D by 4-D
a6  = a(:,:,1,:); b6  = b(:,:,1,:); c6  = c; % 4-D by 4-D
a7  = a(:,:,:,1); b7  = b(:,:,1,:); c7  = c; % 3-D by 4-D
a8  = a(:,:,1,1); b8  = b(:,:,1,:); c8  = c; % 2-D by 4-D
a9  = a(:,:,:,:); b9  = b(:,:,:,1); c9  = c; % 4-D by 3-D
a10 = a(:,:,1,:); b10 = b(:,:,:,1); c10 = c; % 4-D by 3-D
a11 = a(:,:,:,1); b11 = b(:,:,:,1); c11 = c; % 3-D by 3-D
a12 = a(:,:,1,1); b12 = b(:,:,:,1); c12 = c; % 2-D by 3-D
a13 = a(:,:,:,:); b13 = b(:,:,1,1); c13 = c; % 4-D by 2-D
a14 = a(:,:,1,:); b14 = b(:,:,1,1); c14 = c; % 4-D by 2-D
a15 = a(:,:,:,1); b15 = b(:,:,1,1); c15 = c; % 3-D by 2-D
a16 = a(:,:,1,1); b16 = b(:,:,1,1); c16 = c; % 2-D by 2-D
for d1 = 1:5
    for d2 = 1:2 
         % The size of C is adjusted by CHECKALLSHIFTS
         c1(:,:,d1,d2) =  a1(:,:,d1,d2) *  b1(:,:,d1,d2);
         c2(:,:,d1,d2) =  a2(:,:, 1,d2) *  b2(:,:,d1,d2);
         c3(:,:,d1,d2) =  a3(:,:,d1, 1) *  b3(:,:,d1,d2);
         c4(:,:,d1,d2) =  a4(:,:, 1, 1) *  b4(:,:,d1,d2);
         c5(:,:,d1,d2) =  a5(:,:,d1,d2) *  b5(:,:, 1,d2);
         c6(:,:,d1,d2) =  a6(:,:, 1,d2) *  b6(:,:, 1,d2);
         c7(:,:,d1,d2) =  a7(:,:,d1, 1) *  b7(:,:, 1,d2);
         c8(:,:,d1,d2) =  a8(:,:, 1, 1) *  b8(:,:, 1,d2);
         c9(:,:,d1,d2) =  a9(:,:,d1,d2) *  b9(:,:,d1, 1);
        c10(:,:,d1,d2) = a10(:,:, 1,d2) * b10(:,:,d1, 1);
        c11(:,:,d1,d2) = a11(:,:,d1, 1) * b11(:,:,d1, 1);
        c12(:,:,d1,d2) = a12(:,:, 1, 1) * b12(:,:,d1, 1);
        c13(:,:,d1,d2) = a13(:,:,d1,d2) * b13(:,:, 1, 1);
        c14(:,:,d1,d2) = a14(:,:, 1,d2) * b14(:,:, 1, 1);
        c15(:,:,d1,d2) = a15(:,:,d1, 1) * b15(:,:, 1, 1);
        c16(:,:,d1,d2) = a16(:,:, 1, 1) * b16(:,:, 1, 1);
    end
end

checkallshifts(a1, b1, c1);
checkallshifts(a2, b2, c2);
checkallshifts(a3, b3, c3);
checkallshifts(a4, b4, c4);
checkallshifts(a5, b5, c5);
checkallshifts(a6, b6, c6);
checkallshifts(a7, b7, c7);
checkallshifts(a8, b8, c8);
checkallshifts(a9, b9, c9);
checkallshifts(a10,b10,c10);
checkallshifts(a11,b11,c11);
checkallshifts(a12,b12,c12);
checkallshifts(a13,b13,c13);
checkallshifts(a14,b14,c14);
checkallshifts(a15,b15,c15);
checkallshifts(a16,b16,c16);


%--------------------------------------------------------------------------
function checkallshifts(a,b,c)
% A, B, C are 4-D arrays

id = [1 2]; % Internal dimensions of A, B, C
if size(a,3)==1 && size(b,3)==1, c = c(:,:,1,:); end
if size(a,4)==1 && size(b,4)==1, c = c(:,:,:,1); end
order = [5 4 3 id]; % Turns 4-D array into 5-D (with one leading singleton)
checkallshifts2(a,b,c, order,    [4 5]);
checkallshifts2(a,b,c, [3 id 4], [2 3]);
checkallshifts2(a,b,c, [id 3 4], [1 2]);


%--------------------------------------------------------------------------
function checkallshifts2(a,b,c, order, newIDs)
% A, B, C are 4-D arrays

d1 = newIDs(1);
idA0 = newIDs;
idB0 = newIDs;
a0 = permute(a, order);
b0 = permute(b, order);
c0 = permute(c, order);
sizeA0 = size(a0); 
sizeB0 = size(b0);
onesA0 = find(sizeA0~=1, 1) - 1; % Number of leading singleton dims
onesB0 = find(sizeB0~=1, 1) - 1; % Number of leading singleton dims
onesA0 = min(onesA0, d1-1); % Max allowed value is D1-1 (D1 and D2 are IDs)
onesB0 = min(onesB0, d1-1);

for i = 0:onesA0,      a = shiftdim(a0, i); idA = idA0 - i;
    for j = 0:onesB0,  b = shiftdim(b0, j); idB = idB0 - j;
                       c = shiftdim(c0, min(i,j));
                       checkallsyntaxes(a,b,c, idA, idB);
    end
end


%--------------------------------------------------------------------------
function checkallsyntaxes(a,b,c, idA, idB)

global numberofchecks  min_of_all_e  max_of_all_e  syntaxes

idC = max(idA, idB);
idA1 = idA(1); idA2 = idA(2);
idB1 = idB(1); idB2 = idB(2);
idC1 = idC(1); idC2 = idC(2);
disp ' '
disp '------------------------------------------------------------------'
disp ' '
idstrA(10+idA*6) = '•';
idstrB(10+idB*6) = '•';
idstrC(10+idC*6) = '•';
diffA = idA2 - ndims(a);
diffB = idB2 - ndims(b);
diffC = idC2 - ndims(c);
disp (idstrA), fprintf ('Size of A:'), disp ([size(a) ones(1,diffA)])
disp (idstrB), fprintf ('Size of B:'), disp ([size(b) ones(1,diffB)])
disp (idstrC), fprintf ('Size of C:'), disp ([size(c) ones(1,diffC)])

p = size(a, idA1);
q = size(a, idA2);
r = size(b, idB1);
s = size(b, idB2);
scalarsinA = false;
scalarsinB = false;
if isequal([p q], [1 1]), scalarsinA = true; end
if isequal([r s], [1 1]), scalarsinB = true; end

c1s   = c;
c1xs  = c;
c2a   = c;
c2b   = c;
c3a   = c;
c3b   = c;
c4ab  = c;
c4abs = c;
c4b   = c;
c4c   = c;
Nofsyntaxes = 12;
checklog(1:Nofsyntaxes) = false;

% Syntax 1 (Matrices by Matrices)
checklog(1) = true;    
c1 = multiprod(a,b, idA, idB);  % Mat by Mat long syntax
if isequal(idA, idB)
    checklog(2) = true;    
    c1s  = multiprod(a,b, idA); % Mat by Mat short syntax
    if idA1==1
        checklog(3) = true;    
        c1xs  = multiprod(a,b); % Mat by Mat extra-short syntax
    end
end
    
if p == 1, Ad2 = remove1(a, idA1); end % 1 x Q
if q == 1, Ad1 = remove1(a, idA2); end % P x 1
if r == 1, Bd2 = remove1(b, idB1); end % 1 x S
if s == 1, Bd1 = remove1(b, idB2); end % R x 1

% REMOVING ONLY OUTER SINGLETONS IS ALWAYS POSSIBLE
% REMOVING INNER SINGLETONS IS NOT ALWAYS POSSIBLE
% NOTE: either BOTH or NONE inner dimensions can be singletons    

% Syntax 2 (removing right outer singleton)
if s == 1
    if scalarsinB % Syntax 2b (multiprod returns 2-D)
        checklog(5) = true;    
                     c2b = multiprod(a, Bd1, [idA1 idA2], idB1);
    else % Syntax 2a
        checklog(4) = true;    
        c2a = insert1(multiprod(a, Bd1, [idA1 idA2], idB1), idC2);
    end
end

% Syntax 3 (removing left outer singleton)
if p == 1 
    if scalarsinA % Syntax 3b (multiprod returns 2-D)
        checklog(7) = true;    
                     c3b = multiprod(Ad2, b, idA1, [idB1 idB2]);
    else % Syntax 3a
        checklog(6) = true;    
        c3a = insert1(multiprod(Ad2, b, idA1, [idB1 idB2]), idC1);
    end
end

% Syntaxes 4a, 4b (removing BOTH outer singletons)
if all([p     s] == 1)
    % Syntax 4a/b
    if scalarsinA && scalarsinB, dimaddC = idC1; % In both    - C = 1x1
    elseif scalarsinA,           dimaddC = idC2; % Only in A  - C = Rx1
    elseif scalarsinB,           dimaddC = idC1; % Only in B  - C = 1xQ
    else                         dimaddC = idC1; % No scalars - C = 1x1
    end
    checklog(8) = true;    
    c4ab = insert1(multiprod(Ad2, Bd1, idA1, idB1), dimaddC);
    if idA1 == idB1
        checklog(9) = true;
        c4abs = insert1(multiprod(Ad2, Bd1, idA1), dimaddC);
    end
    
    % Syntax 4b - Inner product
    %             C is 1x1
    if q == r
        checklog(10) = true;    
        c4b = insert1(multiprod(Ad2, Bd1, [0 idA1], [idB1 0]), idC1);
    end

end

% Syntax 4c - Outer product (removing BOTH inner singletons)
%             C is PxS
if all([  q r  ] == 1)
    checklog(11) = true;    
    c4c = multiprod(Ad1, Bd2, [idA1 0], [0 idB1]); 
end

% Checking results
checkresult = [ isequal(c1, c), ...
                isequal(c1s, c), ...
                isequal(c1xs, c), ...
                isequal(c2a, c) , ...
                isequal(c2b, c) , ...
                isequal(c3a, c) , ...
                isequal(c3b, c) , ...
                isequal(c4ab, c), ...
                isequal(c4abs, c), ...
                isequal(c4b, c), ...
                isequal(c4c, c)];
            
% Displaying results
temp = (1 : Nofsyntaxes) * 6;
logstring(temp(checklog)) = '*';
disp ' ', disp 'The checked syntaxes are indicated by ''*'':'
disp (syntaxes)
disp (logstring)

mine = 0;
maxe = 0;
if  all(checkresult)
    disp ' ', disp 'The multiproduct is correct :-)'
else
    dotstring(temp(~checkresult)) = '•';    
    disp (dotstring)
    disp ' ', disp 'Warning:'
    disp '    The syntaxes indicated by ''•'' returned wrong multiproducts :-('
    disp '    (probably just a rounding error; see min and max error below)'
    e1 = c1 - c;
    e2 = c1s  - c;
    e3 = c1xs - c;
    e4 = c2a - c;
    e5 = c2b - c;
    e6 = c3a - c;
    e7 = c3b - c;
    e8 = c4ab - c;
    e9 = c4abs - c;
    e10 = c4b - c;
    e11 = c4c - c;
    e_vector = [e1(:); e2(:); e3(:); e4(:); e5(:); e6(:); e7(:); ...
                e8(:); e9(:); e10(:); e11(:)];
    mine = min(e_vector);
    maxe = max(e_vector);
    disp ' '
    disp ( ['Minimun error:    ' num2str(mine)] )
    disp ( ['Maximum error:    ' num2str(maxe)] )
    disp ( ['MATLAB precision: ' num2str(eps)] )
%     disp ' '
%     disp 'Press any key to continue...'
%     pause
end

min_of_all_e = min(mine, min_of_all_e);
max_of_all_e = max(maxe, max_of_all_e);
numberofchecks(checklog) = numberofchecks(checklog) + 1;



%--------------------------------------------------------------------------
function b = insert1(a, newdim)
%Adding a singleton dimension to an array.

if newdim <= ndims(a)
    sizA = size(a);
    sizB = [sizA(1:newdim-1), 1, sizA(newdim:end)];
    b = reshape(a, sizB);
else 
    b = a;
end



%--------------------------------------------------------------------------
function b = remove1(a, dim)
%Removing a singleton dimension from an array.

if dim < ndims(a)
    sizA = size(a);
    sizB = [sizA(1:dim-1), sizA(dim+1:end), 1];
    b = reshape(a, sizB);
else 
    b = a;
end