% Demonstrates additiv plotting
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2019-02-28
%
% Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.
%

function []=gui()
    addpath('ffmatlib');

    mesh_file = 'demo_Lshape_mesh.msh';
    fe_space_file = 'demo_Lshape_vh.txt';
    data_file = 'demo_Lshape_data.txt';
    handles.view = [-166,40];

%     mesh_file = 'capacitor_2d.msh';
%     fe_space_file = 'capacitor_vh_2d.txt';
%     data_file = 'capacitor_data_2d.txt';
%     handles.view = [-47,24];

    [data.p,data.b,data.t,data.nv,data.nbe,data.nt,data.labels]=ffreadmesh(mesh_file);
    [data.vh]=ffreaddata(fe_space_file);
    [data.u,data.Ex,data.Ey]=ffreaddata(data_file);
    handles.data = data;
    winWidth=700;
    winHeight=470;
    handles.f=figure('Name','ffpdeplot() demo','NumberTitle','off','MenuBar','none','visible', 'off', ...
                     'Position', [75 75 winWidth winHeight],'Resize','off');
    set(handles.f,'color',[0.95,0.95,0.95]);
    uicontrol('style','text','String', 'Chose Configuration:','position', [30 (60+30*12) 150 20], ...
              'HorizontalAlignment','left','BackgroundColor',[0.95,0.95,0.95],'FontWeight','bold','FontSize',9);
    handles.names = {'CenteredColor','Patch','Patch+Mesh','Contour','Quiver','Mesh','Boundary', ...
                     '3D-Patch','3D-Patch+Mesh','3D-Mesh','Grid'};
    for k = 1:numel(handles.names)
        str = handles.names{k};
        handles.cb(k) = uicontrol('Style','checkbox','Position',[30 (55+30*k) 130 20], ...
                                  'String',str,'Callback',@checkbox_callback);
    end
    set(handles.cb, 'BackgroundColor', [0.95,0.95,0.95]);
    uicontrol('Style', 'pushbutton', 'String','Close', 'Position', ...
                  [30 30 90 23],'Callback', @pushbuttonExit_Callback);
    %Use subplots if colorbars are involved
    handles.subplot=subplot('Position',[0.3 0.1 0.65 0.8]);
    %handles.ax=axes('Units','pixels','Position',[230,50,350,320]);
    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'VhSeq',handles.data.vh, ...
              'XYData',handles.data.u,'Mesh','on', 'Boundary','on','ColorMap',jet(100), ...
              'BDLabels',handles.data.labels,'BDShowText','on');
    axis tight;
    title('Patch+Mesh Boundary');
    set(handles.cb(3),'Value',1);
    set(handles.cb(7),'Value',1);
    figure(handles.f);
    drawnow;
    guidata(handles.f, handles);
    uiwait(handles.f);
    gh=ishghandle(handles.f);
    if gh
        close(handles.f);
    end
end

function checkbox_callback(hObject, eventdata)
    handles = guidata(hObject);
    %cla;
    %delete the subplot to get rid of the colorbar if there is any
    delete(handles.subplot);
    handles.subplot=subplot('Position',[0.3 0.1 0.65 0.8]);
    hold on;
    val=get(handles.cb,'Value');
    titlestr=' ';
    colrange='minmax';
    for i = 1:numel(val)
        if (val{i})
            switch (get(handles.cb(i),'String'))
                case ('CenteredColor')
                    colrange='centered';
                case ('Mesh')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'Mesh','on', ...
                              'Boundary','off');
                    titlestr=sprintf('%s %s',titlestr,'Mesh');
                case ('Patch')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'VhSeq',handles.data.vh, ...
                              'XYData',handles.data.u,'Mesh','off', 'Boundary','off','ColorMap',jet(100), ...
                              'ColorRange',colrange);
                    titlestr=sprintf('%s %s',titlestr,'Patch');
                case ('Patch+Mesh')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'VhSeq',handles.data.vh, ...
                              'XYData',handles.data.u,'Mesh','on', 'Boundary','off','ColorMap',jet(100), ...
                              'ColorRange',colrange);
                    titlestr=sprintf('%s %s',titlestr,'Patch+Mesh');
                case ('Contour')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'VhSeq',handles.data.vh, ...
                              'XYData',handles.data.u, 'Mesh','off', 'Boundary','off', 'Contour','on', ...
                              'CColor','b','XYStyle','off','CGridParam',[150, 150],'ColorBar','off', ...
                              'ColorMap','off','ColorRange',colrange);
                    titlestr=sprintf('%s %s',titlestr,'Contour');
                case ('Quiver')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'VhSeq',handles.data.vh, ...
                              'Mesh','off','Boundary','off','XYStyle','off', ...
                              'ColorBar','off','FlowData',[handles.data.Ex,handles.data.Ey],'FGridParam',[24, 24]);
                    titlestr=sprintf('%s %s',titlestr,'Quiver');
                case ('Boundary')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'Boundary','on', ...
                              'BDLabels',handles.data.labels,'BDShowText','on');
                    titlestr=sprintf('%s %s',titlestr,'Boundary');
                case ('3D-Patch')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'VhSeq',handles.data.vh, ...
                              'XYData',handles.data.u,'ZStyle','continuous','Mesh','off','ColorMap',jet(100), ...
                              'ColorRange',colrange);
                    lighting gouraud;
                    view(handles.view);
                    camlight('headlight');
                    titlestr=sprintf('%s %s',titlestr,'3D-Patch');
                case ('3D-Patch+Mesh')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'VhSeq',handles.data.vh, ...
                              'XYData',handles.data.u,'ZStyle','continuous','Mesh','on','ColorMap',jet(100), ...
                              'ColorRange',colrange);
                    lighting gouraud;
                    view(handles.view);
                    camlight('headlight');
                    titlestr=sprintf('%s %s',titlestr,'3D-Patch+Mesh');
                case ('3D-Mesh')
                    ffpdeplot(handles.data.p,handles.data.b,handles.data.t,'VhSeq',handles.data.vh, ...
                              'XYData',handles.data.u,'ZStyle','continuous','Mesh','on','XYStyle','off', ...
                              'ColorBar','off','ColorMap','off','ColorRange',colrange);
                    view(handles.view);
                    titlestr=sprintf('%s %s',titlestr,'3D-Mesh');
                case ('Grid')
                    grid off;
                    grid on;
                otherwise
                    error('unknown field');
            end
        end
    end
    title(titlestr);
    axis tight;
    ylabel('y');
    xlabel('x');
    guidata(hObject, handles);
end

function pushbuttonExit_Callback(hObject, eventdata)
    handles=guidata(hObject);
    uiresume(handles.f);
end
