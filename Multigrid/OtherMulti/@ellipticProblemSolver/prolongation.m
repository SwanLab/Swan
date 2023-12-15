function d_h = prolongation(d_2h)    
    i_h = @(i_2h) 2*i_2h - 1;
    j_h = @(j_2h) 2*j_2h - 1;
    
    p = [1/4, 1/2, 1/4;
         1/2,   1, 1/2;
         1/4, 1/2, 1/4];
    
    N_2h = size(d_2h, 1);
    
    N_h = i_h(N_2h);
    N_h = j_h(N_2h);    
    
    d_h = zeros(N_h, N_h);
    for j_2h = 1:N_2h
        for i_2h = 1:N_2h
            d_h_loc = d_2h(j_2h, i_2h)*p;
            setNeighbors(j_h(j_2h)-1, i_h(i_2h)-1, d_h_loc(1, 1));                        
            setNeighbors(j_h(j_2h)-1, i_h(i_2h),   d_h_loc(1, 2));
            setNeighbors(j_h(j_2h)-1, i_h(i_2h)+1, d_h_loc(1, 3));
            setNeighbors(j_h(j_2h),   i_h(i_2h)-1, d_h_loc(2, 1)); 
            setNeighbors(j_h(j_2h),   i_h(i_2h),   d_h_loc(2, 2));
            setNeighbors(j_h(j_2h),   i_h(i_2h)+1, d_h_loc(2, 3));
            setNeighbors(j_h(j_2h)+1, i_h(i_2h)-1, d_h_loc(3, 1));
            setNeighbors(j_h(j_2h)+1, i_h(i_2h),   d_h_loc(3, 2));
            setNeighbors(j_h(j_2h)+1, i_h(i_2h)+1, d_h_loc(3, 3));
        end
    end   
    
    function setNeighbors(j, i, val)
        if j >= 1 && j <= N_h && i >= 1 && i <= N_h
            d_h(j, i) = d_h(j, i) + val;
        end
    end
end

