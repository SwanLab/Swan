function fh = openfig_opt( fname )

if exist(fname,'file')
    fh = openfig(fname);
    % mp = get(0, 'MonitorPositions');
    % if size(mp,1) > 1
    %     set(fh,'Position',mp(2,:));
    % end
else
    fh = 'Figure file not found!';
end

end

