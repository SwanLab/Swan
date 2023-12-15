function v = JacobiIter(A, v0, g)
    N = size(v0, 1);
    v = zeros(size(v0));
    v(1, :)   = g(1,   :);
    v(end, :) = g(end, :);
    v(:, 1)   = g(:,   1);
    v(:, end) = g(:, end);
	
	v_prev = v0;
    % indices of inner nodes
    i_in = 2:N-1;
    j_in = 2:N-1;
    v(j_in, i_in) = (g(j_in, i_in) -...
                   A.d.*v_prev(j_in,   i_in-1)-...
                   A.c.*v_prev(j_in-1, i_in)-...
                   A.e.*v_prev(j_in+1, i_in)-...
                   A.b.*v_prev(j_in,   i_in+1))./A.a;
end