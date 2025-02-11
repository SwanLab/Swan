function list = updateList(path)
    init_list = dir(path);
    list = init_list;
    
    k = 0;
    for i = 1:length(list)
        if strcmpi(init_list(i).name,'.') || strcmpi(init_list(i).name,'..')
            list(i-k) = [];
            k = k+1;
        end
    end
end

