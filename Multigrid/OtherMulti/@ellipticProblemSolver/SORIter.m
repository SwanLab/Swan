function v = SORIter(A, v_prev, g, omega)
    N = size(v_prev, 1);
     
    v = zeros(size(v_prev));
    v(1, :)   = g(1,   :);
    v(end, :) = g(end, :);
    v(:, 1)   = g(:,   1);
    v(:, end) = g(:, end);
	
    for j = 2:N-1
        for i = 2:N-1
            v(j, i) = v_prev(j, i) + omega*((g(j, i)-...
              A.d(j-1, i-1)*v(j, i-1)-...
              A.c(j-1, i-1)*v(j-1, i)-...
              A.e(j-1, i-1)*v_prev(j+1, i)-...
              A.b(j-1, i-1)*v_prev(j, i+1))/A.a(j-1, i-1)-v_prev(j, i));
        end
    end

end