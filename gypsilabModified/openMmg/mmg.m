classdef mmg < handle
% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique  
% Algiane Froehly, CARDAMOME, INRIA-SOFT 
% LGPL Lesser General Public License v3.0. 
% Remeshing using Mmg tools : https://www.mmgtools.org             
properties (SetAccess = public, GetAccess = private)
    oldFileName = 'original';
    newFileName = 'refined';    
end

properties 
    % Edge 
    Aniso    = [];   % AUTOMATIC ANISOTROPY
    Optim    = [];   % AUTOMATIC MESH IMPROVEMENT
    Hgrad    = [];   % MAXIMAL RATIO BETEEN 2 ADJACENT EDGES
    Hgradreq = [];   % RATIO BETWEEN REQUIRED ENTITIES AND NEIGHBOURG 
    Hmax     = [];   % MAXIMAL EDGE LENGTH
    Hmin     = [];   % MINIMAL EDGE LENGTH
    Hsiz     = [];   % CONSTANT EDGE LENGTH
    Hausd    = [];   % HAUSSDORF DISTANCE
    
    % Adaptation
    Nofem    = [];   % REMESHING WITHOUT FINITE ELEMENT ENFORCEMENT 
    Noinsert = [];   % CONSTANT NODES NUMBER
    Nomove   = [];   % NO POINT RELOCATION 
    Nosurf   = [];   % NO SURFACE MODIFICATION  
    Noswap   = [];   % NO EDGE SWAPPING
    Angle    = [];   % SHARP ANGLE DETECTION (°) 
    
    % General
    Memory  = [];    % MAXIMUM MEMORY USED (M0)
    Verbose = [];    % VERBOSITY
    Options = [];    % PRINT OPTIONS DEFINITION
    
    % Mesh data
    Mesh    = [];    % MESH OBJECT FROM GYPSILAB (vtx, elt, col)
    Req     = [];    % REQUIRED ENTITIES AT ELEMENTS
    Map     = [];    % SIZE MAP AT VERTICES
    
    % File names
    oldMeshFileName
    newMeshFileName
    levelSetOldFileName
    levelSetNewFileName
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    function obj = mmg(varargin)
        % Initialize
        if (nargin==0)
            
        % Input is mesh    
        elseif (nargin==1)
            obj.Mesh = varargin{1};
            
        % Input is mesh and Haussdorf distance
        elseif (nargin==2)
            obj.Mesh  = varargin{1};
            obj.Hausd = ['-hausd ',num2str(varargin{2})];
            
        else
            error('mmg.m : unavailable case')
        end
                
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function aniso(mmg)
        mmg.Aniso = '-A';
    end
    
    function optim(mmg)
        mmg.Optim = '-optim';
    end
    
    function hgrad(mmg,val)
        mmg.Hgrad = ['-hgrad ',num2str(val)];
    end
    
    function hgradreq(mmg,val)
        mmg.Hgradreq = ['-hgradreq ',num2str(val)];
    end
    
    function hmax(mmg,val)
        mmg.Hmax = ['-hmax ',num2str(val)];
    end
    
    function hmin(mmg,val)
        mmg.Hmin = ['-hmin ',num2str(val)];
    end
    
    function hsiz(mmg,val)
        mmg.Hsiz  = ['-hsiz ',num2str(val)];
    end
    
    function hausd(mmg,val)
        mmg.Hausd  = ['-hausd ',num2str(val)];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADAPTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function nofem(mmg)
        mmg.Nofem = '-nofem';
    end
    
    function noinsert(mmg)
        mmg.Noinsert = '-noinsert';
    end
    
    function nomove(mmg)
        mmg.Nomove = '-nomove';
    end
    
    function nosurf(mmg)
        mmg.Nosurf = '-nosurf';
    end
    
    function noswap(mmg)
        mmg.Noswap = '-noswap';
    end
    
    function angle(mmg,val)
        mmg.Angle = ['-ar ',num2str(val)];
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function memory(mmg,val)
        mmg.Memory = ['-m ',num2str(val)];
    end
    
    function verbose(mmg,val)
        mmg.Verbose = ['-v ',num2str(val)];
    end
    
    function options(mmg)
        mmg.Options = '-h ';
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MESH DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function mesh(mmg,val)
        mmg.Mesh = val;
    end
    
    function req(mmg,val)
        mmg.Req = val;
    end
    
    function map(mmg,val)
        mmg.Map = val;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function mesh = run(mmg)
        mesh = mmgRun(mmg);
    end
    
    function mesh = runLs(mmg)
        mesh = mmgRunLs(mmg);
    end    
    
    function v = get.oldMeshFileName(obj)
      v = [obj.oldFileName,'.mesh'];
    end
    
    function v = get.newMeshFileName(obj)
       v = [obj.newFileName,'.mesh'];
    end
    
    function v = get.levelSetOldFileName(obj)
        v = [obj.oldFileName,'.sol'];
    end
    
    function v = get.levelSetNewFileName(obj)
        v = [obj.newFileName,'.sol'];
    end
    
end
end
