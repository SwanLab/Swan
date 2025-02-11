function message = sysrequirements_for_testing(varargin)
% Checking whether the functions needed for testing are installed.
%     In order to run m-functions timing_arraylab_engines and 
%     testing_memory_usage, you need to have installed in your system:
%         - MATLAB R2007a or later (you need the builtin function BSXFUN)
%         - GENOP (MATLAB Central file #10333)
%         - BSXFUN substitute (MATLAB Central file #23005)
%         - TIMEIT (MATLAB Central file #18798)

% Paolo de Leva
% University of Rome, Foro Italico, Rome, Italy
% 2009 Feb 22

disp ' '
disp(['You are running MATLAB version ' version])
errorflag = false;

  bsxfunflag = false;
  bsxmexflag = false;
   genopflag = false;
bsxtimesflag = false;
  timeitflag = false;

for i = 1:nargin
    switch varargin{i}
        case 'bsxfun',     bsxfunflag = true;
        case 'bsxmex',     bsxmexflag = true;
        case 'genop',       genopflag = true;
        case 'bsxtimes', bsxtimesflag = true;
        case 'timeit',     timeitflag = true;
    end
end

% Check to see if BSXFUN is a builtin function.
if bsxfunflag && ~exist('bsxfun', 'builtin')
    disp ' '
    w = (['WARNING: You need MATLAB 7.4 (R2007a) or a later version, because  '
          '         you need the builtin function BSXFUN (its replacements for'
          '         releases prior to R2007a are also tested here).           ']);
    disp (w)
    errorflag = true;
end
            
if genopflag && ~exist('genop', 'file')
    disp ' '
    w = (['WARNING: You need to download GENOP (MATLAB Central file #10333,'
          '         http://www.mathworks.com/matlabcentral/fileexchange)   '
          '         and move it into the current directory.                ']);
    disp (w)
    errorflag = true;
end

if bsxmexflag && (~exist('bsxfun', 'file') || ~exist('bsx_times', 'file'))
    disp ' '
    w = (['WARNING: You need to download Schwarz''s BSXFUN substitute (file #23005'
          '         http://www.mathworks.com/matlabcentral/fileexchange), install'
          '         it, add to the MATLAB search path the directory in which you '
          '         installed it (use PATHTOOL to set the path).                 '
          '         Do not forget to run MAKE_BSX_MEX to build BSX_TIMES.        ']);
      disp (w)
    errorflag = true;
end

if bsxtimesflag && ~exist('bsx_times', 'file')
    disp ' '
    w = (['WARNING: You need to download Schwarz''s BSXFUN substitute (file #23005'
          '         http://www.mathworks.com/matlabcentral/fileexchange), install'
          '         it, delete BSXFUN.M to avoid a conflict with the builtin     '
          '         BSXFUN, run MAKE_BSX_MEX to build BSX_TIMES, and move        '
          '         BSX_TIMES to the current directory.                          ']);
      disp (w)
    errorflag = true;
end

if timeitflag
    verdate = datenum(version('-date'), 'mmmm dd, yyyy');
    minimumdate = datenum(2005,6,21); % Correct version of functions TIC/TOC
    if verdate < minimumdate
        disp ' '
        w = (['WARNING: You need MATLAB 7.1 (R14SP3) or a later version, because '
              '         you need an accurate and precise version of the builtin  '
              '         functions TIC and TOC, supporting the undocumented syntax'
              '         T = TIC; ...; ELAPSED_TIME = TOC(T).                     ']);
        disp (w)
        errorflag = true;
    elseif ~exist('timeit', 'file')
        disp ' '
        w = (['WARNING: You need to download TIMEIT (MATLAB Central file #18798,   '
              '         http://www.mathworks.com/matlabcentral/fileexchange),      '
              '         install it, and add to the MATLAB search path the directory'
              '         in which you installed it (use PATHTOOL to set the path).  ']);
        disp (w)
        errorflag = true;
    end
end

if errorflag
    message = (['Before running this function you must install the above\n', ... 
                'mentioned software']);
else
    message = false;
end


