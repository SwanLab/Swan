function adjust_figure_properties(hfig, font_size, picturewidth, hw_ratio)

%   ADJUST_FIGURE_PROPERTIES Sets properties for the figure appearance and export settings.
%
%   Inputs:     hfig        - Handle to the figure (this specifies which figure to adjust).
%               font_size   - Desired font size for all labels and axes.
%               picturewidth - Width of the figure in centimeters.
%               hw_ratio    - Height-to-width ratio to maintain figure proportions.
%
%   Description: This function standardizes figure properties by setting font sizes for labels, enabling LaTeX
%               interpretation for text, and adjusting the figure’s dimensions. It also prepares the figure
%               for export.

    % Set font sizes and other properties:

    set(findall(hfig, '-property', 'FontSize'), 'FontSize', font_size);                        % Adjust font size
    set(findall(hfig, '-property', 'Box'), 'Box', 'on');                                       % Set box around the plot
    set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex');                    % Use LaTeX interpreter
    set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');  % Use LaTeX for tick labels
    
    legend_handle = findobj(hfig, 'Type', 'Legend');   % Find legend objects in the figure
    if ~isempty(legend_handle)
        set(legend_handle, 'Box', 'on');              % Turn off the legend box
    end

    set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio * picturewidth]);  % Set figure size
    pos = get(hfig, 'Position');  % Get current position
    set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)]);  % Set paper size
end