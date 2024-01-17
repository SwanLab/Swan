% Loops through recursive V-cycles while the error norm is within 10e-7
function [y_new, vcycle] = MG_full_cycle(nx_current, ny_current, A, f)
r = ones(size(f));
y_old = zeros(size(f));
vcycle = 0;
while(norm(r) > 10.0e-7)
    y_new =  MG(nx_current, ny_current, A, f, y_old);
    y_old = y_new;
    r = A*y_new - f;
    vcycle = vcycle + 1;
end
end