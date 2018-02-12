function close_fig_except(figname)

% Closes all figures except the one specified in figname

%% Get figure handles
fh = findall(0,'type','figure'); % gets all figure handles

%% Get names
nfig = length(fh);
name = cell(nfig,1);

for i = 1:nfig
    name{i} = get(fh(i),'name');
end

%% Filter by name
idx = ~ismember(name,figname);
fh = fh(idx);

%% Close figures
close(fh);

end
