function make_weights_evolution_gif(historyW, qLATENT, VOL, LOCAL, historyTitle)
% Create either:
%   (A) a GIF of the weights' evolution (default), or
%   (B) an interactive FIG with a slider to scrub frames (if LOCAL.INTERACTIVE = true)
if nargin == 0
    load('tmp2.mat')
    LOCAL.INTERACTIVE = true;
    LOCAL.MAKE_GIF = false;
    
else
    if nargin < 5 || isempty(historyTitle), historyTitle = {}; end
if nargin < 4 || isempty(LOCAL), LOCAL = struct; end
if ~isfield(LOCAL,'SHOW_LEGEND'), LOCAL.SHOW_LEGEND = true; end
if ~isfield(LOCAL,'MAKE_GIF'),     LOCAL.MAKE_GIF = true; end
if ~isfield(LOCAL,'GIF_FILE'),     LOCAL.GIF_FILE = 'weights_evolution.gif'; end
if ~isfield(LOCAL,'DELAY_TIME'),   LOCAL.DELAY_TIME = 0.25; end
if ~isfield(LOCAL,'DPI'),          LOCAL.DPI = 120; end
if ~isfield(LOCAL,'INTERACTIVE'),  LOCAL.INTERACTIVE = true; end

    
end

if ~isfield(LOCAL,'SAVE_FIG'),  LOCAL.SAVE_FIG = false; end
if ~isfield(LOCAL,'FIG_FILE'),  LOCAL.FIG_FILE = 'weights_evolution.fig'; end

[nPoints, nClusters] = size(historyW{1});
nFrames = numel(historyW);

% Figure / axes
vis = 'off'; if LOCAL.INTERACTIVE, vis = 'on'; end
f = figure('Color','w','Visible',vis);
ax = axes('Parent',f); hold(ax,'on'); box(ax,'on');
xlabel(ax,'qLATENT'); ylabel(ax,'weights (non-zero) / VOL * 100');

% Lines
hLines = gobjects(nPoints,1);
for iii = 1:nPoints
    hLines(iii) = plot(ax, qLATENT, nan(1,nClusters), ...
        'DisplayName', sprintf('w_{%d}',iii));
end
if LOCAL.SHOW_LEGEND, legend(ax,'show','Location','eastoutside'); end

% Y-limits from all frames
maxVal = 0;
for it = 1:nFrames
    maxVal = max(maxVal, max(historyW{it}(:)));
end
ylim(ax, [0, max(1e-12, maxVal/VOL*100)*1.05]);

% Store shared state (for callbacks)
S.historyW = historyW; S.qLATENT = qLATENT; S.VOL = VOL; S.LOCAL = LOCAL;
S.nPoints = nPoints; S.nClusters = nClusters; S.nFrames = nFrames;
S.f = f; S.ax = ax; S.hLines = hLines; S.historyTitle = historyTitle;
guidata(f,S);

if LOCAL.INTERACTIVE
    % Layout + slider UI
    f.Position(3:4) = [1000 650];
    ax.Position = [0.08 0.20 0.80 0.75];
    
    sld = uicontrol('Style','slider','Units','normalized',...
        'Min',1,'Max',max(1,nFrames),'Value',1,...
        'SliderStep', stepFromFrames(nFrames),...
        'Position',[0.08 0.10 0.80 0.05],...
        'Callback',@onSlide);
    
    txt = uicontrol('Style','text','Units','normalized',...
        'Position',[0.90 0.10 0.08 0.05],...
        'HorizontalAlignment','left','String','1');
    
    S.sld = sld; S.txt = txt; guidata(f,S);
    updateFrame(f,1);
    set(f,'Visible','on');
    if LOCAL.SAVE_FIG
        try
            savefig(f, LOCAL.FIG_FILE);
            fprintf('FIG written to: %s\n', LOCAL.FIG_FILE);
        catch ME
            warning('Could not save FIG: %s', ME.message);
        end
    end
    return; % no GIF in interactive mode
end

% --- GIF mode ---
f.Position(3:4) = [1000 600];
for it = 1:nFrames
    updateFrame(f,it);
    frame = getframe(f);
    [imind, cm] = rgb2ind(frame2im(frame), 256);
    if it == 1
        imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', 'DelayTime', LOCAL.DELAY_TIME);
    else
        imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', 'WriteMode','append', ...
            'DelayTime', LOCAL.DELAY_TIME);
    end
end
close(f);
fprintf('GIF written to: %s\n', LOCAL.GIF_FILE);

% ===== Nested helpers (keep inside this function to avoid nesting errors) ====

    function updateFrame(fig, it)
        S = guidata(fig);
        W = S.historyW{it};
        for jj = 1:S.nPoints
            y = W(jj,:)/S.VOL*100;
            if all(W(jj,:)==0)
                set(S.hLines(jj), 'YData', nan(1,S.nClusters), ...
                    'Visible','off','HandleVisibility','off');
            else
                set(S.hLines(jj), 'YData', y, ...
                    'Visible','on','HandleVisibility','on', ...
                    'DisplayName', sprintf('w_{%d}', jj));
            end
        end
        if S.LOCAL.SHOW_LEGEND, legend(S.ax,'show'); end
        if it <= numel(S.historyTitle) && ~isempty(S.historyTitle{it})
            title(S.ax, S.historyTitle{it});
        else
            title(S.ax, sprintf('Weights vs qLATENT — frame %d/%d', it, S.nFrames));
        end
        drawnow;
    end

    function onSlide(src,~)
        fig = ancestor(src,'figure');
        S = guidata(fig);
        val = max(1, min(S.nFrames, round(get(src,'Value'))));
        set(src,'Value',val);
        if isfield(S,'txt') && isgraphics(S.txt)
            set(S.txt, 'String', sprintf('%d',val));
        end
        updateFrame(fig,val);
    end

    function step = stepFromFrames(nF)
        if nF <= 1
            step = [1 1];
        else
            step = [1/(nF-1), min(1, 5/(nF-1))];
        end
    end
end




% function make_weights_evolution_gif(historyW, qLATENT, VOL, LOCAL, historyTitle)
% % Create a GIF that shows the evolution of weights per integration point
% % historyW: cell, each cell is [nPoints x nClusters] weights at an iteration
% % qLATENT: 1 x nClusters
% % VOL: scalar to normalize (% of volume)
% % LOCAL: settings (GIF_FILE, DELAY_TIME, DPI, SHOW_LEGEND)
% % historyTitle: cell of strings (titles per frame)
%
% [nPoints, nClusters] = size(historyW{1});
%
% % Prepare figure
% f = figure('Color','w','Visible','off');
% ax = axes('Parent',f); hold(ax,'on'); box(ax,'on');
% xlabel(ax,'qLATENT'); ylabel(ax,'weights (non-zero) / VOL * 100');
%
% % Pre-create line handles for consistency across frames
% hLines = gobjects(nPoints,1);
% for iii = 1:nPoints
%     hLines(iii) = plot(ax, qLATENT, nan(1,nClusters), 'DisplayName', sprintf('w_{%d}',iii));
% end
% if LOCAL.SHOW_LEGEND
%     legend(ax,'show','Location','eastoutside');
% end
%
% % axis limits based on all frames to avoid resizing
% maxVal = 0;
% for it = 1:numel(historyW)
%     maxVal = max(maxVal, max(historyW{it}(:)));
% end
% ylim(ax, [0, max(1e-12, maxVal/VOL*100)*1.05]);
%
% % Write GIF
% for it = 1:numel(historyW)
%     W = historyW{it};
%     for iii = 1:nPoints
%         y = W(iii,:)/VOL*100;
%         if all(W(iii,:)==0)
%             set(hLines(iii), 'YData', nan(1,nClusters), 'Visible','off', 'HandleVisibility','off');
%         else
%             set(hLines(iii), 'YData', y, 'Visible','on', 'HandleVisibility','on', ...
%                 'DisplayName', sprintf('w_{%d}', iii));
%         end
%     end
%     if LOCAL.SHOW_LEGEND
%         legend(ax,'show');
%     end
%     if it <= numel(historyTitle)
%         title(ax, historyTitle{it});
%     else
%         title(ax, sprintf('Weights vs qLATENT — frame %d', it));
%     end
%
%     drawnow;
%
%     % Capture frame and write/append to GIF
%     f.Position(3:4) = [1000 600]; % make it reasonably wide
%     frame = getframe(f);
%     [imind, cm] = rgb2ind(frame2im(frame), 256);
%     if it == 1
%         %             imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', ...
%         %                 'LoopCount', Inf, 'DelayTime', LOCAL.DELAY_TIME);
%         %     imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', ...
%         % 'LoopCount', 1, 'DelayTime', LOCAL.DELAY_TIME);
%
%         imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', 'DelayTime', LOCAL.DELAY_TIME);
%     else
%         imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', ...
%             'WriteMode', 'append', 'DelayTime', LOCAL.DELAY_TIME);
%     end
% end
% close(f);
% fprintf('GIF written to: %s\n', LOCAL.GIF_FILE);
% end