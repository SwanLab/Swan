function v = SeidelIter(A, v0, g)
    N = size(v0, 1);
   
    v = zeros(size(v0));
    v(1, :)   = g(1,   :);
    v(end, :) = g(end, :);
    v(:, 1)   = g(:,   1);
    v(:, end) = g(:, end);
	
	v_prev = v0;    
    for j = 2:N-1
        for i = 2:N-1                    
            v(j, i) = (g(j, i)-...
              A.d(j-1, i-1)*v(j, i-1)-...
              A.c(j-1, i-1)*v(j-1, i)-...
              A.e(j-1, i-1)*v_prev(j+1, i)-...
              A.b(j-1, i-1)*v_prev(j, i+1))/A.a(j-1, i-1);
        end
    end   
end