%% Creates contour(), quiver() as well as patch() plots from FreeFem++ 2D simulation data
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-06-15
%
% [handles,varargout] = ffpdeplot (points, boundary, triangles, varargin)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  ffpdeplot() is a function specially tailored to FreeFem++ that offers most
%  of the features of the classic Matlab pdeplot() command. contour() plots
%  (2D iso values), quiver() plots (2D vector fields) and patch() plots
%  (2D map data) can be created as well as their combinations. In addition
%  domain border edges can be selectively displayed and superimposed to the plot
%  data. The FEM mesh is entered through its vertices, the boundary values and
%  the triangles as provided by the FreeFem++ command
%  savemesh(Th, "filename.msh"). The finite element connectivity data as well
%  as the PDE simulation data are provided using the FreeFem++ macros
%  ffExportVh(filename.txt, Th, Vh) and ffExportData1(filename.txt, u).
%

%% Input Parameters
%
%   [handles,varargout] = ffpdeplot (p,b,t,'PARAM1',val1,'PARAM2',val2, ...)
%
%   specifies parameter name/value pairs to control the input file format
%
%      Parameter       Value
%      'VhSeq'        Finite element connectivity
%                        FreeFem++ macro definition
%      'XYData'       PDE data used to create the plot
%                        FreeFem++ macro definition
%      'XYStyle'      Coloring choice
%                        'interp' (default) | 'off'
%      'ZStyle'       Draws 3D surface plot instead of flat 2D Map plot
%                        'continuous' | 'off' (default)
%      'ColorMap'     ColorMap value or matrix of such values
%                        'off' | 'cool' (default) | colormap name | three-column matrix of RGB triplets
%      'ColorBar'     Indicator in order to include a colorbar
%                        'on' (default) | 'off' | 'northoutside' ...
%      'CBTitle'      Colorbar Title
%                        (default=[])
%      'ColorRange'   Range of values to adjust the colormap thresholds
%                        'off' | 'minmax' (default) | 'centered' | 'cropminmax' | 'cropcentered' | [min,max]
%      'Mesh'         Switches the mesh off / on
%                        'on' | 'off' (default)
%      'MColor'       Color to colorize the mesh
%                        'auto' (default) | RGB triplet | 'r' | 'g' | 'b'
%      'RLabels'      Meshplot of specified regions
%                        [] (default) | [region1,region2,...]
%      'RColors'      Colorize regions with a specific color (linked to 'RLabels')
%                        'b' (default) | three-column matrix of RGB triplets
%      'Boundary'     Shows the domain boundary / edges
%                        'on' | 'off' (default)
%      'BDLabels'     Draws boundary / edges with a specific label
%                        [] (default) | [label1,label2,...]
%      'BDColors'     Colorize boundary / edges with a specific color (linked to 'BDLabels')
%                        'r' (default) | three-column matrix of RGB triplets
%      'BDShowText'   Shows the labelnumber on the boundary / edges
%                        'on' | 'off' (default)
%      'BDTextSize'   Size of labelnumbers on the boundary / edges
%                        scalar value greater than zero
%      'BDTextWeight' Character thickness of labelnumbers on the boundary / edges
%                        'normal' (default) | 'bold'
%      'Contour'      Isovalue plot
%                        'off' (default) | 'on'
%      'CStyle'       Contour plot style
%                        'solid' (default) | 'dashed' | 'dashedneg'
%      'CColor'       Isovalue color (can be monochrome or flat)
%                        'flat' | [0,0,0] (default) | RGB triplet, three-element row vector | 'r' | 'g' | 'b'
%      'CLevels'      Number of isovalues used in the contour plot
%                        (default=10)
%      'CGridParam'   Number of grid points used for the contour plot
%                        'auto' (default) | [N,M]
%      'Title'        Title
%                        (default=[])
%      'XLim'         Range for the x-axis
%                        'minmax' (default) | [min,max]
%      'YLim'         Range for the y-axis
%                        'minmax' (default) | [min,max]
%      'ZLim'         Range for the z-axis
%                        'minmax' (default) | [min,max]
%      'DAspect'      Data unit length of the xy- and z-axes
%                        'off' | 'xyequal' (default) | [ux,uy,uz]
%      'FlowData'     Data for quiver plot
%                        FreeFem++ point data | FreeFem++ triangle data
%      'FColor'       Color to colorize the quiver arrows
%                        'b' (default) | RGB triplet | 'r' | 'g'
%      'FGridParam'   Number of grid points used for quiver plot
%                        'auto' (default) | [N,M]

%% Output Parameters
%
%  List of handles to the graphic objects.
%

%% Licence
%
%  Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
%  This program is free software: you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see
%  <https://www.gnu.org/licenses/>.
%

%% Code
function [hh,varargout] = ffpdeplot(points,boundary,triangles,varargin)

    if (nargin<3)
        printhelp();
        error('wrong number arguments');
    end

    numvarargs = length(varargin);
    optsnames = {'VhSeq', 'XYData', 'XYStyle', 'Mesh', 'MColor', 'Boundary', 'BDLabels', ...
                 'BDColors', 'BDShowText', 'BDTextWeight', 'BDTextSize', ...
                 'RLabels', 'RColors', 'ZStyle', 'Contour', ...
                 'CGridParam', 'CColor', 'CLevels', 'CStyle', ...
                 'ColorMap', 'ColorBar', 'CBTitle', 'ColorRange', ...
                 'Title', 'XLim', 'YLim', 'ZLim', 'DAspect', ...
                 'FlowData', 'FColor', 'FGridParam'};

    vararginval = {[], [], 'interp', 'off', [], 'off', [], ...
                   'r', 'off', [], [], ...
                   [], 'b', 'off', 'off', ...
                   'auto', [0,0,0], 10, 'solid', ...
                   'cool', 'on', [], 'minmax', ...
                   [], 'minmax', 'minmax', 'minmax', 'xyequal', ...
                   [], 'b', 'auto'};

    if (numvarargs>0)
        if (~mod(numvarargs,2))
            for i=1:2:(numvarargs-1)
                pos=find(strcmpi(optsnames,varargin(i)));
                if ~isempty(pos)
                    vararginval(pos)=varargin(i+1);
                else
                    printhelp();
                    fprintf('%s\n',char(varargin(i)));
                    error('unknown input parameter');
                end
            end
        else
            printhelp();
            error('wrong number arguments');
        end
    end

    [vhseq, xyrawdata, xystyle, showmesh, mcolor, showboundary, bdlabels, ...
     bdcolors, bdshowtext, bdtxtweight, bdtxtsize, ...
     rlabels, rcolors, zstyle, contourplt, ...
     cgridparam, ccolor, isolevels, contourstyle, ...
     setcolormap, showcolbar, colorbartitle, colorrange, ...
     plottitle, plotxlim, plotylim, plotzlim, axisaspect, ...
     flowrawdata, fcolor, fgridparam] = vararginval{:};

    %newplot() checks the values of the NextPlot-properties and prepare the
    %figure for plotting based on these values. If there is no current figure,
    %newplot creates one.
    hax=newplot();
    fig=get(hax,'Parent');
    oldnextplotval{1}=get(hax,'nextplot');
    oldnextplotval{2}=get(fig,'nextplot');
    %switch hold on
    set(hax,'nextplot','add');
    set(fig,'nextplot','add');
    %fprintf('nextplot hax: %s, fig: %s\n',oldnextplotval{1},oldnextplotval{2});

    %Used for various output handles: ColorBar, Patch, Contour, Quiver
    hh=[];
    %Used for Contour labels, if there are any
    varargout{1}=[];
    varargout{2}=[];
    is2dmode=true;
    %points=rowvec(points);
    %triangles=rowvec(triangles);
    [xmesh,xpts,ymesh,ypts] = prepare_mesh(points,triangles);
    if ~isempty(xyrawdata)
        xyrawdata=rowvec(xyrawdata);
        %Decode the rawdata based on the Finite Element Space
        [elementType,xydata]=convert_pde_data(points,triangles,vhseq,xyrawdata);
        %Based on the Element Space create various Meshes for plotting and interpolation
        [xdata,xdataz,ydata,ydataz]=prepare_coordinates(elementType,xmesh,ymesh,true);
        %Convert the color data according corresponding to xdata, xdataz and interpolation
        [cdata,cdataz,cdatainterpn]=prepare_data(elementType,triangles,xydata,xmesh,ymesh,true);
        %3D and 2D Plots invoking patch() command
        if strcmpi(contourplt,'off')
            %2D Map/Density Plot
            if strcmpi(zstyle,'off')
                %A colored Plot
                if (strcmpi(xystyle,'interp'))
                    %With Mesh
                    if ~strcmpi(showmesh,'off')
                        meshcol = getmcol(mcolor,[0 0 0]);
                        switch (elementType)
                            case {'P0','P1'}
                                %No subtriangles - mesh and data is equal
                                hh=patch(xdata,ydata,cdata,'EdgeColor',meshcol);
                            case {'P1b','P2'}
                                %Plot the data - subtriangles (refined mesh)
                                hh=patch(xdata,ydata,cdata,'EdgeColor','none');
                                %Plot the mesh - defined by triangle vertices
                                %The plots must ne split otherwise parts of
                                %the mesh can be obscured by the subtriangledata
                                hhm=patch(xmesh,ymesh,[1 1 1],'FaceColor','none', ...
                                          'EdgeColor',meshcol);
                                hh=[hh; hhm];
                            otherwise
                                error('Unknown Lagrangian Finite Element. Only P1, P1b and P2 allowed');
                        end
                    else
                        %Without Mesh
                        hh=patch(xdata,ydata,cdata,'EdgeColor','none');
                    end
                %Uncolored Plot (Mesh only).
                %Note: Case is duplicate. In this case coloring is explicetly
                %disabled by the XYStyle parameter.
                else
                    meshcol = getmcol(mcolor,[0 0 1]);
                    hh=patch(xmesh,ymesh,[1 1 1],'FaceColor','none', ...
                             'EdgeColor',meshcol);
                end
                %Once xyrawdata is given there is no way to create an empty plot
                view(2);
            %3D Surface Plots
            else
                is2dmode=false;
                %A colored Map / Surface Plot
                if (strcmpi(xystyle,'interp'))
                    %With Mesh
                    if ~strcmpi(showmesh,'off')
                        meshcol = getmcol(mcolor,[0 0 0]);
                        switch (elementType)
                            case {'P0','P1'}
                                hh=patch(xdata,ydata,cdata,cdata,'EdgeColor',meshcol);
                            case {'P1b','P2'}
                                hh=patch(xdata,ydata,cdata,cdata,'EdgeColor','none');
                                %For P2: Display the mesh at the triangle points {1,6,2,4,3,5}
                                %otherwise parts of mesh can be obscured by the previous patch plot data
                                %For P1b: Plot the triangle vertices = mesh
                                hhm=patch(xdataz,ydataz,cdataz,[1 1 1],'FaceColor','none','EdgeColor',meshcol);
                                hh=[hh; hhm];
                            otherwise
                                error('Unknown Lagrangian Finite Element. Only P1, P1b and P2 allowed');
                        end
                    else
                        %Without Mesh
                        hh=patch(xdata,ydata,cdata,cdata,'EdgeColor','none');
                    end
                %Uncolored 3D Plot - Mesh only
                else
                    meshcol = getmcol(mcolor,[0 0 1]);
                    switch (elementType)
                        case {'P0','P1'}
                            hh=patch(xdata,ydata,cdata,[1 1 1],'FaceColor','none', 'EdgeColor',meshcol);
                        case {'P1b','P2'}
                            hh=patch(xdataz,ydataz,cdataz,[1 1 1],'FaceColor','none', 'EdgeColor',meshcol);
                        otherwise
                            error('Unknown Lagrangian Finite Element. Only P1, P1b and P2 allowed');
                    end
                end
                %Once xyrawdata is given there is no way to create an empty plot
                if (~strcmpi(plotzlim,'minmax'))
                     zlim(plotzlim);
                     zf=1.3*(max(plotzlim)-min(plotzlim));
                else
                     zf=1.3*(max(max(cdata))-min(min(cdata)));
                end
                if strcmpi(axisaspect,'xyequal')
                    yf=(max(max(ydata))-min(min(ydata)));
                    xf=(max(max(xdata))-min(min(xdata)));
                    daspect([max(xf,yf) max(xf,yf) zf]);
                else
                    if isnumeric(axisaspect)
                        daspect(axisaspect);
                    end
                end
                view(3);
            end
        %Contour Plot
        else
            if (~strcmpi(cgridparam,'auto') && isnumeric(cgridparam))
                N=cgridparam(1);
                M=cgridparam(2);
            else
                [~,nt]=size(triangles);
                N=sqrt(nt);
                M=N;
            end
            if strcmpi(plotylim,'minmax')
                ymin = min(min(ydata));
                ymax = max(max(ydata));
            else
                ymin = min(plotylim);
                ymax = max(plotylim);
            end
            if strcmpi(plotxlim,'minmax')
                xmin = min(min(xdata));
                xmax = max(max(xdata));
            else
                xmin = min(plotxlim);
                xmax = max(plotxlim);
            end

            x=linspace(xmin,xmax,N);
            y=linspace(ymin,ymax,M);
            [X,Y]=meshgrid(x,y);

            if exist('fftri2gridfast','file')
                switch (elementType)
                    case {'P0','P1'}
                        C=fftri2gridfast(X,Y,xdata,ydata,cdata);
                    case {'P1b','P2'}
                        C=fftri2gridfast(X,Y,xmesh,ymesh,cdatainterpn);
                    otherwise
                        error('Unknown Lagrangian Finite Element. Only P0, P1, P1b and P2 allowed');
                end
            else
                if (N>100)
                    fprintf('Note: To improve runtime build MEX function fftri2gridfast() from fftri2gridfast.c\n');
                end
                switch (elementType)
                    case {'P0','P1'}
                        C=fftri2grid(X,Y,xdata,ydata,cdata);
                    case {'P1b','P2'}
                        C=fftri2grid(X,Y,xmesh,ymesh,cdatainterpn);
                    otherwise
                        error('Unknown Lagrangian Finite Element. Only P0, P1, P1b and P2 allowed');
                end
            end

            if (strcmpi(xystyle,'interp'))
                hh=patch(xdata,ydata,cdata,'EdgeColor','none');
            else
                hh=[];
            end
            switch (contourstyle)
                %Dashed style iso value plot
                case ('dashed')
                    %Workaround Matlab R2013
                    if strcmpi(ccolor,'flat')
                        [clab,hhc]=contour(X,Y,C,isolevels,'--');
                    else
                        [clab,hhc]=contour(X,Y,C,isolevels,'--','LineColor',ccolor);
                    end
                    varargout{1}=clab;
                %Same as before but positive values are displayed with solid lines
                %and negative values are displayed with dashed lines
                case ('dashedneg')
                    if (length(isolevels)==1)
                        step = (max(max(C))-min(min(C)))/(isolevels+1);
                        isolevels = min(min(C))+step:step:max(max(C))-step;
                    end
                    isopos = isolevels(isolevels>=0);
                    if (length(isopos)==1)
                        isopos = [isopos isopos];
                    end;
                    isoneg = isolevels(isolevels<0);
                    if(length(isoneg)==1)
                        isoneg = [isoneg isoneg];
                    end
                    %Workaround Matlab R2013
                    if strcmpi(ccolor,'flat')
                        [clabneg,hhcneg]=contour(X,Y,C,isoneg);
                        [clabpos,hhcpos]=contour(X,Y,C,isopos);
                    else
                        [clabneg,hhcneg]=contour(X,Y,C,isoneg,'--','LineColor',ccolor);
                        [clabpos,hhcpos]=contour(X,Y,C,isopos,'LineColor',ccolor);
                    end
                    hhc=[hhcneg;hhcpos];
                    varargout{1}=clabneg;
                    varargout{2}=clabpos;
                %Default is a solid line plot
                otherwise
                    %Workaround Matlab R2013
                    if strcmpi(ccolor,'flat')
                        [clab,hhc]=contour(X,Y,C,isolevels);
                    else
                        [clab,hhc]=contour(X,Y,C,isolevels,'LineColor',ccolor);
                    end
                    hh=[hh; hhc];
                    varargout{1}=clab;
            end
            %Add patch plot handle if there is patch plot data
            if ~isempty(hh)
                hh=[hh; hhc];
            else
                hh=hhc;
            end
            view(2);
        end
        %set colormap only for the axes we are currently plotting to and not for
        %the whole figure. does not work for older matlab versions
        if ~(strcmpi(setcolormap,'off'))
            colormap(hax,setcolormap);
        end
        if ~(strcmpi(colorrange,'off'))
            if (isnumeric(colorrange))
                if (min(colorrange) == max(colorrange))
                   error('''ColorRange'': Must be a numeric 2-element vector where LIM1 < LIM2');
                end
                caxis(colorrange);
            else
                if strcmpi(colorrange,'minmax') || strcmpi(colorrange,'centered')
                    %Set to [min max] of the whole mesh
                    if strcmpi(colorrange,'minmax')
                        caxis([min(min(cdata)) max(max(cdata))]);
                    else
                    %Colormap symmetric around zero ('centered')
                        caxis([-max(max(abs(cdata))) max(max(abs(cdata)))]);
                    end
                else
                    %Set to [min max] of cropped area (auto ranging)
                    [~,sz2]=size(cdata);
                    keep = true(1,sz2);
                    if ~strcmpi(plotxlim,'minmax')
                        keepx = (xdata > min(plotxlim)) & (xdata < max(plotxlim));
                        keep = keep & (keepx(1,:) & keepx(2,:) & keepx(3,:));
                    end
                    if ~strcmpi(plotylim,'minmax')
                        keepy = (ydata > min(plotylim)) & (ydata < max(plotylim));
                        keep = keep & (keepy(1,:) & keepy(2,:) & keepy(3,:));
                    end
                    crng = [min(min(cdata(:,keep))) max(max(cdata(:,keep)))];
                    if (isempty(crng) || (min(crng) == max(crng)))
                        %E.g. too much magnification
                        error('''ColorRange'': No color spread in the cropped section / try ''minmax''');
                    end
                    if strcmpi(colorrange,'cropminmax')
                        caxis(crng);
                    else
                        %Colormap symmetric around zero ('cropcentered')
                        caxis([-max(abs(crng)) max(abs(crng))]);
                    end
                end
            end
        end
        if ~(strcmpi(showcolbar,'off'))
            if strcmpi(showcolbar,'on')
                hcb=colorbar;
            else
                hcb=colorbar(showcolbar);
            end
            hh=[hcb; hh];
            if ~isempty(colorbartitle)
                title(hcb,colorbartitle);
            end
        end
    %Uncolored, flat 2D mesh plot (no xyrawdata is given)
    else
        if ~strcmpi(showmesh,'off')
            meshcol = getmcol(mcolor,[0 0 1]);
            if ~isempty(rlabels)
                if isnumeric(rcolors)
                    %bdcolors as numeric RGB data (three column RGB triplets)
                    [nRColors,nRGB] = size(rcolors);
                    if (nRGB ~= 3)
                        error('''RColors'': Must be a three column matrix of RGB triplets');
                    end
                else
                    %as colorstring, e.g. 'r' or 'red'
                    nRColors = 1;
                end
                for i=1:numel(rlabels)
                    tmp=triangles(1:3,(triangles(4,:)==rlabels(i)));
                    if any(tmp)
                        xreg=[xpts(tmp(1,:)); xpts(tmp(2,:)); xpts(tmp(3,:))];
                        yreg=[ypts(tmp(1,:)); ypts(tmp(2,:)); ypts(tmp(3,:))];
                        if (nRColors > 1)
                            hq=patch(xreg,yreg,[1 1 1],'FaceColor','none','EdgeColor',rcolors(i,:));
                        else
                            hq=patch(xreg,yreg,[1 1 1],'FaceColor','none','EdgeColor',meshcol);
                        end
                        if ~isempty(hh)
                            hh=[hh; hq];
                        else
                            hh=hq;
                        end
                    else
                        fprintf('Region Label:%i\n',rlabels(i));
                        error('Region Label not found in Mesh');
                    end
                end
            else
                hh=patch(xmesh,ymesh,[1 1 1],'FaceColor','none','EdgeColor',meshcol);
            end
            view(2);
        end
        %If "mesh" is switched of we can exit at this point without any plot
    end

    %Quiver
    if ~isempty(flowrawdata)
        flowrawdata=rowvec(flowrawdata);
        if (~isempty(xyrawdata) && ~(isequal(size(xyrawdata), size(flowrawdata(1,:)))))
            error('''XYData'' has different size than ''FlowData'': To plot different size data use incremental plotting.');
        end
        [elementType,flowdata]=convert_pde_data(points,triangles,vhseq,flowrawdata);
        [xdata,~,ydata,~]=prepare_coordinates(elementType,xmesh,ymesh,false);
        [udata,vdata,~,~,udatainterpn,vdatainterpn]=prepare_data(elementType,triangles,flowdata,xmesh,ymesh,false);
        if (~strcmpi(fgridparam,'auto') && isnumeric(fgridparam))
            N=fgridparam(1);
            M=fgridparam(2);
        else
            N=20;
            M=20;
        end
        if strcmpi(plotylim,'minmax')
            ymin = min(min(ydata));
            ymax = max(max(ydata));
        else
            ymin = min(plotylim);
            ymax = max(plotylim);
        end
        if strcmpi(plotxlim,'minmax')
            xmin = min(min(xdata));
            xmax = max(max(xdata));
        else
            xmin = min(plotxlim);
            xmax = max(plotxlim);
        end

        x=linspace(xmin,xmax,N);
        y=linspace(ymin,ymax,M);
        [X,Y]=meshgrid(x,y);

        if exist('fftri2gridfast','file')
            switch (elementType)
                case {'P0','P1'}
                    [U,V]=fftri2gridfast(X,Y,xdata,ydata,udata,vdata);
                case {'P1b','P2'}
                    [U,V]=fftri2gridfast(X,Y,xmesh,ymesh,udatainterpn,vdatainterpn);
                otherwise
                    error('Unknown Lagrangian Finite Element. Only P0, P1, P1b and P2 allowed');
            end
        else
            switch (elementType)
                case {'P0','P1'}
                    [U,V]=fftri2grid(X,Y,xdata,ydata,udata,vdata);
                case {'P1b','P2'}
                    [U,V]=fftri2grid(X,Y,xmesh,ymesh,udatainterpn,vdatainterpn);
                otherwise
                    error('Unknown Lagrangian Finite Element. Only P0, P1, P1b and P2 allowed');
            end
        end

        idx=(~isnan(U)) & (~isnan(V));
        hq=quiver(X(idx),Y(idx),U(idx),V(idx),'Color',fcolor);
        if ~isempty(hh)
            hh=[hh; hq];
        else
            hh=hq;
        end
    end

    %Adds the domain boundary (border) to all plot types, or creates a border plot
    %if nothing has been drawn yet
    if ~strcmpi(showboundary,'off')
        boundary=rowvec(boundary);
        %Check for the number of boundary colors and do first checks
        if isnumeric(bdcolors)
            %bdcolors as numeric RGB data (three column RGB triplets)
            [nBDColors,nRGB] = size(bdcolors);
            if (nRGB ~= 3)
                error('''BDColors'': Must be a three column matrix of RGB triplets');
            end
        else
            %bdcolors as colorstring, e.g. 'r' or 'red'
            nBDColors = 1;
        end
        if ~isempty(bdlabels)
            %Check if there are same number of bdcolors than labels
            %If only one color value is given then color everything by this color
            if (numel(bdlabels) ~= nBDColors) && (nBDColors > 1)
                error('''BDColors'': Number of BDColors must be equal to number of BDLabels');
            end
            for i=1:numel(bdlabels)
                keep=(boundary(3,:)==bdlabels(i));
                if any(keep)
                    if (nBDColors > 1)
                        line([xpts(boundary(1,keep));xpts(boundary(2,keep))], ...
                             [ypts(boundary(1,keep));ypts(boundary(2,keep))],'Color',bdcolors(i,:),'LineWidth',2);
                    else
                        line([xpts(boundary(1,keep));xpts(boundary(2,keep))], ...
                             [ypts(boundary(1,keep));ypts(boundary(2,keep))],'Color',bdcolors,'LineWidth',2);
                    end
                    if strcmpi(bdshowtext,'on')
                        textpos=find(keep,1,'first');
                        txt=text(xpts(boundary(1,textpos)),ypts(boundary(1,textpos)),num2str(bdlabels(i)));
                        if ~isempty(bdtxtweight)
                            set(txt,'FontWeight',bdtxtweight);
                        end
                        if ~isempty(bdtxtsize)
                            set(txt,'FontSize',bdtxtsize);
                        end
                    end
                else
                    fprintf('Boundary Label:%i\n',bdlabels(i));
                    error('Boundary Label not found in Mesh');
                end
            end
        else
            %No BDLabels specified: All labels are displayed without text
            if (nBDColors == 1)
                line([xpts(boundary(1,:));xpts(boundary(2,:))], ...
                     [ypts(boundary(1,:));ypts(boundary(2,:))],'Color',bdcolors,'LineWidth',2);
            else
                error('''BDColors'': Multple BDColors only with BDLabel specifier possible');
            end
        end
    end

    if ~isempty(plottitle)
        title(plottitle);
    end
    if ~strcmpi(plotxlim,'minmax')
        xlim(plotxlim);
    end
    if ~strcmpi(plotylim,'minmax')
        ylim(plotylim);
    end
    if (is2dmode) && (~strcmpi(axisaspect,'off'))
        if strcmpi(axisaspect,'xyequal')
            daspect([1 1 1]);
        else
            if isnumeric(axisaspect)
                daspect(axisaspect);
            end
        end
    end
    %Restore axes and figure handle properties
    set(hax,'nextplot',oldnextplotval{1});
    set(fig,'nextplot',oldnextplotval{2});
end

function [S] = rowvec(S)
    [sz1,sz2]=size(S);
    if sz1>sz2
        S=S.';
    end
end

function meshcol = getmcol(mcolor,defaultcol)
    if ~isempty(mcolor) && ~strcmpi(mcolor,'auto')
        meshcol=mcolor;
    else
        meshcol=defaultcol;
    end
end

%Display helpscreen
function printhelp()
    fprintf('%s\n\n','Invalid call to ffpdeplot. Correct usage is:');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (p,b,t,varargin)');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (p,b,t,''Boundary'',''on'')');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (p,b,t,''Boundary'',''on'',''Mesh'',''on'')');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (p,b,t,''VhSeq'',vh,''XYData'',u)');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (p,b,t,''VhSeq'',vh,''XYData'',u,''Contour'',''on'')');
    fprintf('%s\n',' -- [handles,varargout] = ffpdeplot (p,b,t,''VhSeq'',vh,''FlowData'',v,''Boundary'',''on'')');
    fprintf('\n');
    fprintf('''VhSeq''        Finite element connectivity\n');
    fprintf('''XYData''       PDE data used to create the plot\n');
    fprintf('''XYStyle''      Coloring choice (default=''interp'')\n');
    fprintf('''ZStyle''       Draws 3D surface plot instead of flat 2D Map plot (default=''off'')\n');
    fprintf('''ColorMap''     ColorMap value or matrix of such values (default=''on'')\n');
    fprintf('''ColorBar''     Indicator in order to include a colorbar\n');
    fprintf('''CBTitle''      Colorbar Title (default=[])\n');
    fprintf('''ColorRange''   Range of values to adjust the color thresholds (default=''minmax'')\n');
    fprintf('''Mesh''         Switches the mesh off / on (default=''off'')\n');
    fprintf('''MColor''       Color to colorize the mesh (default=''auto'')\n');
    fprintf('''RColors''      Colorize regions with a specific color (linked to ''RLabels'')\n');
    fprintf('''RLabels''      Meshplot of specified regions\n');
    fprintf('''Boundary''     Shows the domain boundary / edges (default=''off'')\n');
    fprintf('''BDLabels''     Draws boundary / edges with a specific label (default=[])\n');
    fprintf('''BDColors''     Colorize boundary / edges with color (default=''r'')\n');
    fprintf('''BDShowText''   Shows the labelnumber on the boundary / edges (default=''off'')\n');
    fprintf('''BDTextSize''   Size of labelnumbers on the boundary / edges\n');
    fprintf('''BDTextWeight'' Character thickness of labelnumbers on the boundary / edges\n');
    fprintf('''Contour''      Isovalue plot (default=''off'')\n');
    fprintf('''CStyle''       Contour line style (default=''patch'')\n');
    fprintf('''CColor''       Isovalue color (default=[0,0,0])\n');
    fprintf('''CLevels''      Number of isovalues used in the contour plot (default=10)\n');
    fprintf('''CGridParam''   Number of grid points used for the contour plot (default=''off'')\n');
    fprintf('''Title''        Title (default=[])\n');
    fprintf('''XLim''         Range for the x-axis (default=''minmax'')\n');
    fprintf('''YLim''         Range for the y-axis (default=''minmax'')\n');
    fprintf('''ZLim''         Range for the z-axis (default=''minmax'')\n');
    fprintf('''DAspect''      Data unit length of the xy- and z-axes (default=''xyequal'')\n');
    fprintf('''FlowData''     Data for quiver plot\n');
    fprintf('''FColor''       Color to colorize the quiver arrows (default=''b'')\n');
    fprintf('''FGridParam''   Number of grid points used for quiver plot (default=''off'')\n');
    fprintf('\n');
end
