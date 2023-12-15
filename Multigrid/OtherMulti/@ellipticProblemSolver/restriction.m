function d_2h = restriction(d_h)
    r = [1/16, 1/8, 1/16;
          1/8, 1/4,  1/8;
         1/16, 1/8, 1/16];

    N_h = size(d_h, 1);
    
    i_h = 1:2:N_h;
    j_h = 1:2:N_h;
    N_2h = length(j_h);
    
    d_2h = zeros(N_2h); 
    d_2h(:, 1) = d_h(j_h, 1);
    d_2h(:, N_2h) = d_h(j_h, N_h);
    d_2h(1, :) = d_h(1, i_h);
    d_2h(N_2h, :) = d_h(N_h, i_h);
    for j_2h = 2:N_2h-1
        for i_2h = 2:N_2h-1
            d_h_loc = getNeighbors(j_h(j_2h),i_h(i_2h));
            d_2h(j_2h, i_2h) = sum(sum(r.*d_h_loc));
        end
    end

    function d_h_loc = getNeighbors(j,i)
        d_h_loc = [d_h(j-1, i-1), d_h(j-1, i), d_h(j-1, i+1);
                   d_h(j,   i-1), d_h(j,   i), d_h(j,   i+1);
                   d_h(j+1, i-1), d_h(j+1, i), d_h(j+1, i+1)];
    end    
end

