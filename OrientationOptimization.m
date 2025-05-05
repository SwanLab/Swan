syms theta real
assume(theta, 'real');

% Define trigonometric terms
c2 = cos(2*theta);
s2 = sin(2*theta);

% Define Rsig matrix
Rsig = [ (1 + c2)/2, (1 - c2)/2,   s2;
         (1 - c2)/2, (1 + c2)/2,  -s2;
         -s2/2,       s2/2,       c2 ];

% Define Reps matrix
Reps = [ (1 + c2)/2, (1 - c2)/2,   s2/2;
         (1 - c2)/2, (1 + c2)/2,  -s2/2;
         -s2,         s2,         c2 ];

% Define a general 3x3 symbolic constitutive matrix C
C = sym('C', [3 3]);
assume(C, 'real'); % optional if you want only real values

% Perform the transformation
C_rotated = simplify(Rsig * C * Reps);

% Display the result
disp('Transformed constitutive matrix C_rotated:');
disp(C_rotated);

%% Example 2
syms theta eps1 eps2 eps12 real
assume(theta, 'real');

% Define trigonometric terms
c2 = cos(2*theta);
s2 = sin(2*theta);

% Define Rsig matrix
Rsig = [ (1 + c2)/2, (1 - c2)/2,   s2;
         (1 - c2)/2, (1 + c2)/2,  -s2;
         -s2/2,       s2/2,       c2 ];

% Define Reps matrix
Reps = [ (1 + c2)/2, (1 - c2)/2,   s2/2;
         (1 - c2)/2, (1 + c2)/2,  -s2/2;
         -s2,         s2,         c2 ];

% Define a general 3x3 symbolic constitutive matrix C
C = sym('C', [3 3]);
assume(C, 'real');  % Optional

% Define strain vector in rotated frame
eps_vec = [eps1; eps2; eps12];

% Compute the rotated form: eps' * Rsig * C * Reps * eps
expr = simplify(eps_vec.' * Rsig * C * Reps * eps_vec);

% Derivative with respect to theta
d_expr_dtheta = simplify(diff(expr, theta));

% Display results
disp('Scalar expression eps'' * Rsig * C * Reps * eps:');
disp(expr);
disp('Derivative with respect to theta:');
disp(d_expr_dtheta);

%% Example 3 
syms theta eps1 eps2 eps12 real
assume(theta, 'real');

% Trig shorthand
c2 = cos(2*theta);
s2 = sin(2*theta);

% Define Rsig matrix
Rsig = [ (1 + c2)/2, (1 - c2)/2,   s2;
         (1 - c2)/2, (1 + c2)/2,  -s2;
         -s2/2,       s2/2,       c2 ];

% Define Reps matrix
Reps = [ (1 + c2)/2, (1 - c2)/2,   s2/2;
         (1 - c2)/2, (1 + c2)/2,  -s2/2;
         -s2,         s2,         c2 ];

% General symmetric 3x3 symbolic constitutive matrix
C = sym('C', [3 3]);
C = (C + C.') / 2;  % Ensure symmetry
assume(C, 'real');

% Strain vector
eps_vec = [eps1; eps2; eps12];

% Scalar expression: eps^T * Rsig * C * Reps * eps
expr = simplify(eps_vec.' * Rsig * C * Reps * eps_vec);

% Derivative with respect to theta
d_expr_dtheta = simplify(diff(expr, theta));

% Substitute to extract coefficients
syms s2 c2 s4 c4 real
subs_expr = subs(d_expr_dtheta, {
    sin(2*theta), cos(2*theta), ...
    sin(4*theta), cos(4*theta)
}, {
    s2, c2, s4, c4
});

% Extract and display coefficients
coeff_s2 = simplify(collect(subs_expr, s2));
coeff_c2 = simplify(collect(subs_expr, c2));
coeff_s4 = simplify(collect(subs_expr, s4));
coeff_c4 = simplify(collect(subs_expr, c4));

disp('Full derivative expression:');
disp(d_expr_dtheta);

disp('--- Coefficients ---');
disp('Coefficient of sin(2θ):');
disp(coeffs(coeff_s2, s2));

disp('Coefficient of cos(2θ):');
disp(coeffs(coeff_c2, c2));

disp('Coefficient of sin(4θ):');
disp(coeffs(coeff_s4, s4));

disp('Coefficient of cos(4θ):');
disp(coeffs(coeff_c4, c4));
