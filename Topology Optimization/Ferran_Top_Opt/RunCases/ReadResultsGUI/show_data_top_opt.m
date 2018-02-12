function show_data_top_opt (handles,results)

% Get selected output option
contents_output = cellstr(get(handles.postprocess_data_menu,'String'));
selected_output = contents_output{get(handles.postprocess_data_menu,'Value')};

%% Visualize UI elements
switch selected_output
    case 'Normal'
        normal_elements(handles,'on');
        word_elements(handles,'off');

    case 'Word'
        normal_elements(handles,'off');
        word_elements(handles,'on');
end

%% Update final results
if exist('results','var')
    numcase = get(handles.case_menu,'Value') - 1;
    set(handles.globaliter_gid,'String',results.niter);
    switch selected_output
        case 'Normal'
            % Update final results
            set(handles.compliance,'String',results.Comp);
            set(handles.complianceP0,'String',results.CompP0);
            set(handles.complianceSIMPALL,'String',results.CompSIMPALL);
            set(handles.complianceThreshold,'String',results.CompThreshold);
            set(handles.V,'String',results.Vol);
            set(handles.P,'String',results.Per);
            set(handles.theta,'String',results.theta);
            set(handles.time,'String',handles.data.RunTimes(numcase));
            set(handles.iter,'String',handles.data.Iterations(numcase));
            set(handles.iterglobal,'String',results.niter);


        case 'Word'
            str = {''};
            % Compliance
            if get(handles.chk_compliance,'Value')
                str = [str;{sprintf('Comp. (real): %.3f',results.Comp)}];
            end
            
            if get(handles.chk_complianceP0,'Value')
                if strcmp(results.CompP0,'-')
                    results.CompP0 = results.Comp;
                end
                str = [str;{sprintf('Compliance: %.3f',results.CompP0)}];
            end
            
            if get(handles.chk_complianceSIMPALL,'Value')
                if ~strcmp(results.CompSIMPALL,'-') % only for SIMP
                    str = [str;{sprintf('Compliance (SIMP-ALL): %.3f',results.CompSIMPALL)}];
                end
            end
            
            if get(handles.chk_complianceThreshold,'Value')
                str = [str;{sprintf('Comp. Threshold: %.3f',results.CompThreshold)}];
            end
            
            % Iterations
            if get(handles.chk_iter,'Value')
                str = [str;{sprintf('Function evaluations: %.0f',handles.data.Iterations(numcase))}];
            end
            
            if get(handles.chk_iterglobal,'Value')
                str = [str;{sprintf('Global iterations: %.0f',results.niter)}];
            end
            
            % Constraints
            if get(handles.chk_V,'Value')
                str = [str;{sprintf('Volume: %.3f',results.Vol)}];
            end
            if get(handles.chk_P,'Value')
                str = [str;{sprintf('Perimeter: %.3f',results.Per)}];
            end
            
            % Performance
            if get(handles.chk_theta,'Value')
                if results.kkttol == 1 % level set
                    str = [str;{sprintf('Theta: %.2f',results.theta)}];
                else % otherwise
                    str = [str;{sprintf('KKT tolerance: %.3e',results.kkttol)}];
                end
            end
            
            if get(handles.chk_time,'Value')
                if iscell(handles.data.RunTimes(numcase))
                    str = [str;{sprintf('Time: %.0fs',handles.data.RunTimes{numcase})}];
                else
                    str = [str;{sprintf('Time: %.0fs',handles.data.RunTimes(numcase))}];
                end
            end
            
            
            
            % Set selected results to text editor
            set(handles.txt_word,'String',str);
            
    end
end


end

function normal_elements(handles,state)

% Modify visible state normal
set(handles.compliance,'Visible',state);
set(handles.complianceP0,'Visible',state);
set(handles.complianceSIMPALL,'Visible',state);
set(handles.complianceThreshold,'Visible',state);
set(handles.V,'Visible',state);
set(handles.P,'Visible',state);
set(handles.theta,'Visible',state);
set(handles.time,'Visible',state);
set(handles.iter,'Visible',state);
set(handles.iterglobal,'Visible',state);

set(handles.txt_compliance,'Visible',state);
set(handles.txt_complianceP0,'Visible',state);
set(handles.txt_complianceSIMPALL,'Visible',state);
set(handles.txt_complianceThreshold,'Visible',state);
set(handles.txt_V,'Visible',state);
set(handles.txt_P,'Visible',state);
set(handles.txt_theta,'Visible',state);
set(handles.txt_time,'Visible',state);
set(handles.txt_iter,'Visible',state);
set(handles.txt_iterglobal,'Visible',state);

end

function word_elements(handles,state)

% Modify visible state word
set(handles.chk_compliance,'Visible',state);
set(handles.chk_complianceP0,'Visible',state);
set(handles.chk_complianceSIMPALL,'Visible',state);
set(handles.chk_complianceThreshold,'Visible',state);
set(handles.chk_V,'Visible',state);
set(handles.chk_P,'Visible',state);
set(handles.chk_theta,'Visible',state);
set(handles.chk_time,'Visible',state);
set(handles.chk_iter,'Visible',state);
set(handles.chk_iterglobal,'Visible',state);

set(handles.txt_word,'Visible',state);
set(handles.copy_txt,'Visible',state);

end
