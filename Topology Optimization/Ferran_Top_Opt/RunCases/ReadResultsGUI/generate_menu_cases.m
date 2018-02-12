function cases = generate_menu_cases( data )

ncases = length(data.NumCase);
cases = cell(ncases+1,1);
cases{1} = '';
for i = 1:ncases
    method = data.Method{i};
    switch method
        case 'SIMP_ALL'
            method = 'SIMP-ALL';
        case 'SIMP'
            method = 'SIMP';
    end
    
    algorithm = data.Algorithm{i};
    switch algorithm
        case 'Projected_gradient'
            algorithm = 'PG';
        case 'level_set'
            algorithm = 'LS';
    end
    
    kernel = data.Kernel{i};
    switch kernel
        case 'PDE'
            kernel = 'PDE';
        case 'P1_kernel'
            kernel = 'P1';
        case 'P0_kernel'
            kernel = 'P0';
    end
    cases{i+1} = sprintf('Case %s - %s %s %s',data.NumCase{i},method,algorithm,kernel);
end

end

