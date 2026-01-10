@ECHO OFF

rem OutputFile: %2\%1.sal

rem set basename = %1
rem set directory = %2
rem set ProblemDirectory = %3

:init_dir
cd %2

move /Y %1.dat %1.cal

%3\AnyExecutable.exe %1.cal

del %1.dic
del outputfile.dat
del outputfile.bon
move /Y outputfile.res %1.flavia.res
