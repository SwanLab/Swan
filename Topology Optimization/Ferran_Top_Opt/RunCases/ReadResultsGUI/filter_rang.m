function index_new = filter_rang(data,index,proccase)

found_match = find(index==true);
k = regexpi(data(1:end-1),proccase);
fcn2modify_ini = find(~cellfun('isempty', k));

% Find function endings
pattern =  'Description';
k = regexpi(data(1:end-1),pattern);
all_fcn_end = find(~cellfun('isempty', k));

% Generate ranges
fcn2modify_end = zeros(size(fcn2modify_ini)); 
for i = 1:length(fcn2modify_ini)
    aux = find(all_fcn_end > fcn2modify_ini(i),1);
    fcn2modify_end(i) = all_fcn_end(aux);
end

% Filter by ranges
index_new = index;
for i = 1:length(found_match)
    aux = found_match(i);
    ub = aux > fcn2modify_ini;
    lb = aux < fcn2modify_end;
    if ~any(ub.*lb) % if it's not inside any valid range, set to false
        index_new(aux) = false;
    end
end

end