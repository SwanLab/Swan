function results = get_opt_results(fh,method)

% Gets the results from the figure handle specified in fh
    
%% Get all axes handles
ax = get(fh,'children'); % gets all axes handles
tit = get(ax,'title'); % gets title handles from axes

% Get title strings
fun = @(x) get(x,'String');
titles = cellfun(fun,tit,'UniformOutput',false);
titles = strrep(titles,'\','');
titles = strrep(titles,'_','');

% Remove unneeded axes
valid = {'Compliance','Volume','Perimeter','theta','theta / mu','Incre gamma'};
idx = ismember(titles,valid);
ax = ax(idx);
titles = titles(idx);

%% Get data from the axes
naxes = length(ax);
for i = 1:naxes
    [x,y,marker] = get_plot_data (ax(i));
    
    ndata = length(marker);
    for j = 1:ndata
        yval = y{j}(end);
        
        switch titles{i}
            % Compliance
            case 'Compliance'
                switch marker{j}
                    case '*' % using threshold
                        results.CompThreshold = yval;

                    case 'x' % real with P0
                        results.CompP0 = yval;

                    case '+' % using SIMP-ALL interpolation (only if SIMP)
                        results.CompSIMPALL = yval;

                    case 'none' % evolution iterations
                        results.Comp = yval;
                end

            % Volume
            case 'Volume'
                results.Vol = yval;

            % Perimeter
            case 'Perimeter'
                results.Per = yval;

            % Theta
            case {'theta','theta / mu'}
                results.theta = yval;
                
            % iepsilon
            case 'iepsilon'
                results.iepsilon = yval;
                
            % incre gamma / kkttol
            case 'Incre gamma'
                results.kkttol = yval;
                                             
        end   
    end
end
results.niter = max(x{1})+1; % goes from 0 to N -> total N+1

% Exceptions or missing data
if ~isfield(results,'CompSIMPALL')
    results.CompSIMPALL = '-';
end
if ~isfield(results,'CompP0')
    results.CompP0 = '-';
end
if ~isfield(results,'CompThreshold')
    results.CompThreshold = '-';
end
if ~isfield(results,'Comp')
    results.Comp = '-';
end
if ~isfield(results,'Vol')
    results.Vol = '-';
end
if ~isfield(results,'Per')
    results.Per = '-';
end
if ~isfield(results,'theta')
    results.theta = '-';
end

if ~isfield(results,'iepsilon')
    results.iepsilon = '-';
end

if ~isfield(results,'kkttol')
    results.kkttol = '-';
end


end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [x,y,marker] = get_plot_data (axh)

lh = get(axh,'Children');  % get all line objects
ndata = length(lh);

x = cell(ndata,1);
y = cell(ndata,1);
marker = cell(ndata,1);

for i = 1:ndata
    x{i} = get(lh(i),'XData');
    y{i} = get(lh(i),'YData');
    marker{i} = get(lh(i),'marker');
end

end