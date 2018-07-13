function varargout = ResultsOptGUI(varargin)
% RESULTSOPTGUI MATLAB code for ResultsOptGUI.fig
%      RESULTSOPTGUI, by itself, creates a new RESULTSOPTGUI or raises the existing
%      singleton*.
%
%      H = RESULTSOPTGUI returns the handle to a new RESULTSOPTGUI or the handle to
%      the existing singleton*.
%
%      RESULTSOPTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTSOPTGUI.M with the given input arguments.
%
%      RESULTSOPTGUI('Property','Value',...) creates a new RESULTSOPTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ResultsOptGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ResultsOptGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ResultsOptGUI

% Last Modified by GUIDE v2.5 03-Mar-2017 14:34:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ResultsOptGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ResultsOptGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ResultsOptGUI is made visible.
function ResultsOptGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ResultsOptGUI (see VARARGIN)

% Choose default command line output for ResultsOptGUI
handles.output = hObject;

% Center figure
movegui(hObject,'center');

% Mark the check box
set(handles.check_close,'Value',1);
set(handles.check_open,'Value',0);
set(handles.check_image,'Value',1);

% Initialize text
handles.menu_str = 'Select a folder first';
handles.dir_str = 'Select a folder to read';
set(handles.case_menu,'String',handles.menu_str);
set(handles.dir_box,'String',handles.dir_str);

% Initialize axes
set(handles.axesim,'Visible','off');

% Set To Word by default
set(handles.postprocess_data_menu,'Value',2);
show_data_top_opt(handles);
set(handles.txt_word,'Max',2); % enable multiple lines

% Default check boxes
set(handles.chk_complianceP0,'Value',1);
set(handles.chk_complianceSIMPALL,'Value',1);
set(handles.chk_V,'Value',1);
set(handles.chk_P,'Value',1);
% set(handles.chk_time,'Value',1);
set(handles.chk_iter,'Value',1);
set(handles.chk_iterglobal,'Value',1);
set(handles.chk_theta,'Value',1);

% Default names
set(handles.chk_compliance,'String',get(handles.txt_compliance,'String'));
set(handles.chk_complianceP0,'String',get(handles.txt_complianceP0,'String'));
set(handles.chk_complianceSIMPALL,'String',get(handles.txt_complianceSIMPALL,'String'));
set(handles.chk_complianceThreshold,'String',get(handles.txt_complianceThreshold,'String'));
set(handles.chk_V,'String',get(handles.txt_V,'String'));
set(handles.chk_P,'String',get(handles.txt_P,'String'));
set(handles.chk_theta,'String',get(handles.txt_theta,'String'));
set(handles.chk_time,'String',get(handles.txt_time,'String'));
set(handles.chk_iter,'String',get(handles.txt_iter,'String'));
set(handles.chk_iterglobal,'String',get(handles.txt_iterglobal,'String'));

% Select paths
handles.user_paths = get_runner_paths;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ResultsOptGUI wait for user response (see UIRESUME)
% uiwait(handles.GUI);


% --- Outputs from this function are returned to the command line.
function varargout = ResultsOptGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in dir_pb.
function dir_pb_Callback(hObject, eventdata, handles)
% hObject    handle to dir_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get updated handles
handles = guidata(hObject);

% Get file and print in box
default_folder = handles.user_paths.default_folder;
path_name = uigetdir(default_folder);

if ~path_name
    return;
end
set(handles.dir_box,'String',path_name);


% Read data in the folder
[handles.data,message,handles.all_iter_saved] = read_files_opt(path_name);
guidata(hObject, handles); % Update handles structure
if ~istable(handles.data)
    set(handles.dir_box,'String',message);
    return;
end

% Set possible cases in menu
cases = generate_menu_cases(handles.data);
set(handles.case_menu,'String',cases);
set(handles.case_menu,'Value',1)

% Reset all data
set(handles.case_box,'String','');
set(handles.method_box,'String','');
set(handles.algorithm_box,'String','');
set(handles.kernel_box,'String','');
set(handles.V_box,'String','');
set(handles.P_box,'String','');
set(handles.lamV_box,'String','');
set(handles.penV_box,'String','');
set(handles.lamP_box,'String','');
set(handles.penP_box,'String','');

set(handles.compliance,'String','');
set(handles.complianceP0,'String','');
set(handles.complianceSIMPALL,'String','');
set(handles.complianceThreshold,'String','');
set(handles.V,'String','');
set(handles.P,'String','');
set(handles.theta,'String','');
set(handles.time,'String','');
set(handles.iter,'String','');
set(handles.iterglobal,'String','');

set(handles.results_file,'String','');
set(handles.globaliter_gid,'String','');
set(handles.txt_word,'String','');

cla(handles.axesim,'reset');
set(handles.axesim,'Visible','off');

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in case_menu.
function case_menu_Callback(hObject, eventdata, handles)
% hObject    handle to case_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get updated handles
handles = guidata(hObject);

% Get selected option
contents = cellstr(get(handles.case_menu,'String'));
selected_case = contents{get(handles.case_menu,'Value')};

if ~isempty(selected_case)
    numcase = get(handles.case_menu,'Value') - 1;
    
    % Show results file and topology
    set(handles.results_file,'String',handles.data.ResultFile{numcase});
    if get(handles.check_image,'Value');
        ResultFile = get(handles.results_file,'String');
        handles.I = postprocessGIDopt(ResultFile,handles.data.Case{numcase},handles.data.Algorithm{numcase},handles.user_paths);
        guidata(hObject, handles); % Update handles structure
        if ~isempty(handles.I)
            set(handles.axesim,'Visible','on');
            axes(handles.axesim);
            imshow(handles.I);
        end
    end
    
    % Open results figure
    if get(handles.check_close,'Value');
        close_fig_except(get(handles.GUI,'Name'));
    end
    fh = openfig_opt(handles.data.FigureFile{numcase});
        
    % Update numerical data
    set(handles.case_box,'String',handles.data.Case{numcase});
    set(handles.method_box,'String',handles.data.Method{numcase});
    set(handles.algorithm_box,'String',handles.data.Algorithm{numcase});
    set(handles.kernel_box,'String',handles.data.Kernel{numcase});
    set(handles.V_box,'String',num2str(handles.data.Vtarget(numcase)));
    set(handles.P_box,'String',num2str(handles.data.Ptarget(numcase)));
    set(handles.lamV_box,'String',num2str(handles.data.lamV(numcase)));
    set(handles.penV_box,'String',num2str(handles.data.penV(numcase)));
    set(handles.lamP_box,'String',num2str(handles.data.lamP(numcase)));
    set(handles.penP_box,'String',num2str(handles.data.penP(numcase)));
    
    % Get final results
    if ishandle(fh)
        handles.results = get_opt_results(fh,handles.data.Method{numcase});
        guidata(hObject, handles); % Update handles structure
        if ~get(handles.check_open,'Value');
            close_fig_except(get(handles.GUI,'Name'));
        end
        show_data_top_opt (handles,handles.results);
        
    else
        set(handles.compliance,'String',fh);
    end
    
else
    cla(handles.axesim,'reset');
    set(handles.axesim,'Visible','off');
    if get(handles.check_close,'Value');
        close_fig_except(get(handles.GUI,'Name'));
    end
    % Reset all data
    set(handles.case_box,'String','');
    set(handles.method_box,'String','');
    set(handles.algorithm_box,'String','');
    set(handles.kernel_box,'String','');
    set(handles.V_box,'String','');
    set(handles.P_box,'String','');
    set(handles.lamV_box,'String','');
    set(handles.penV_box,'String','');
    set(handles.lamP_box,'String','');
    set(handles.penP_box,'String','');

    set(handles.compliance,'String','');
    set(handles.complianceP0,'String','');
    set(handles.complianceSIMPALL,'String','');
    set(handles.complianceThreshold,'String','');
    set(handles.V,'String','');
    set(handles.P,'String','');
    set(handles.theta,'String','');
    set(handles.time,'String','');
    set(handles.iter,'String','');
    set(handles.iterglobal,'String','');

    set(handles.results_file,'String','');
    set(handles.globaliter_gid,'String','');
    set(handles.txt_word,'String','');

    cla(handles.axesim,'reset');
    set(handles.axesim,'Visible','off');
end

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function case_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to case_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function results_file_Callback(hObject, eventdata, handles)
% hObject    handle to results_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of results_file as text
%        str2double(get(hObject,'String')) returns contents of results_file as a double


% --- Executes during object creation, after setting all properties.
function results_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to results_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_close.
function check_close_Callback(hObject, eventdata, handles)
% hObject    handle to check_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_close


% --- Executes on button press in load_topology_pb.
function load_topology_pb_Callback(hObject, eventdata, handles)
% hObject    handle to load_topology_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get updated handles
handles = guidata(hObject);

% Copy value to clipboard
plot_global_iter = get(handles.globaliter_gid,'String');
results_file = get(handles.results_file,'String');


if ~isempty(plot_global_iter) && ~isempty(results_file) 
    numcase = get(handles.case_menu,'Value') - 1;

    [path,file,ext] = fileparts(results_file);
    if str2double(plot_global_iter)+1 > length(handles.all_iter_saved{numcase})
        plot_global_iter = str2double(plot_global_iter);
    else
        plot_global_iter = str2double(plot_global_iter)+1;
    end
    iter = handles.all_iter_saved{numcase}(plot_global_iter);
    file = regexprep(file,'_\d*.flavia',['_',num2str(iter),'.flavia']);
    results_file = [path,'/',file,ext];
    set(handles.results_file,'String',results_file);
    
    ResultFile = results_file;
    handles.I = postprocessGIDopt(ResultFile,handles.data.Case{numcase},handles.data.Algorithm{numcase},handles.user_paths);
    guidata(hObject, handles); % Update handles structure
    if ~isempty(handles.I)
        set(handles.axesim,'Visible','on');
        axes(handles.axesim);
        imshow(handles.I);
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in postprocess_data_menu.
function postprocess_data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to postprocess_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get updated handles
handles = guidata(hObject);

% Show results
if isfield(handles,'results')
    show_data_top_opt (handles,handles.results);
end

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function postprocess_data_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to postprocess_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_pb.
function save_pb_Callback(hObject, eventdata, handles)
% hObject    handle to save_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get updated handles
handles = guidata(hObject);

% Get updated handles
handles = guidata(hObject);

% Get selected option
contents = cellstr(get(handles.case_menu,'String'));
selected_case = contents{get(handles.case_menu,'Value')};

if ~isempty(selected_case)
    numcase = get(handles.case_menu,'Value') - 1;
    results_file = get(handles.results_file,'String');
    
    % Export STL
    exportSTL(results_file,handles.data.Algorithm{numcase},handles.user_paths);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in check_open.
function check_open_Callback(hObject, eventdata, handles)
% hObject    handle to check_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get updated handles
handles = guidata(hObject);

% Get selected option
numcase = get(handles.case_menu,'Value') - 1;

if numcase
    if get(handles.check_open,'Value');
        openfig_opt(handles.data.FigureFile{numcase});
    else
        close_fig_except(get(handles.GUI,'Name'));
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in copyim_pb.
function copyim_pb_Callback(hObject, eventdata, handles)
% hObject    handle to copyim_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get updated handles
handles = guidata(hObject);

% Copy image to clipboard
imclipboard('copy',handles.I);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in chk_compliance.
function chk_compliance_Callback(hObject, eventdata, handles)
% hObject    handle to chk_compliance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_compliance


% --- Executes on button press in chk_complianceP0.
function chk_complianceP0_Callback(hObject, eventdata, handles)
% hObject    handle to chk_complianceP0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_complianceP0


% --- Executes on button press in chk_complianceSIMPALL.
function chk_complianceSIMPALL_Callback(hObject, eventdata, handles)
% hObject    handle to chk_complianceSIMPALL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_complianceSIMPALL


% --- Executes on button press in chk_complianceThreshold.
function chk_complianceThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to chk_complianceThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_complianceThreshold


% --- Executes on button press in chk_V.
function chk_V_Callback(hObject, eventdata, handles)
% hObject    handle to chk_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_V


% --- Executes on button press in chk_P.
function chk_P_Callback(hObject, eventdata, handles)
% hObject    handle to chk_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_P


% --- Executes on button press in chk_theta.
function chk_theta_Callback(hObject, eventdata, handles)
% hObject    handle to chk_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_theta


% --- Executes on button press in chk_time.
function chk_time_Callback(hObject, eventdata, handles)
% hObject    handle to chk_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_time


% --- Executes on button press in chk_iter.
function chk_iter_Callback(hObject, eventdata, handles)
% hObject    handle to chk_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_iter


% --- Executes on button press in chk_iterglobal.
function chk_iterglobal_Callback(hObject, eventdata, handles)
% hObject    handle to chk_iterglobal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_iterglobal



function txt_word_Callback(hObject, eventdata, handles)
% hObject    handle to txt_word (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_word as text
%        str2double(get(hObject,'String')) returns contents of txt_word as a double


% --- Executes during object creation, after setting all properties.
function txt_word_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_word (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in copy_txt.
function copy_txt_Callback(hObject, eventdata, handles)
% hObject    handle to copy_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get updated handles
handles = guidata(hObject);

% Copy value to clipboard
data = get(handles.txt_word,'String');
clipboard('copy',sprintf('%s\n',data{:}));

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in check_image.
function check_image_Callback(hObject, eventdata, handles)
% hObject    handle to check_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_image



function globaliter_gid_Callback(hObject, eventdata, handles)
% hObject    handle to globaliter_gid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of globaliter_gid as text
%        str2double(get(hObject,'String')) returns contents of globaliter_gid as a double


% --- Executes during object creation, after setting all properties.
function globaliter_gid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to globaliter_gid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
