REM @ECHO OFF
REM Identification for arguments
REM basename                          = %1
REM Project directory                 = %2
REM Problem directory                 = %3
 
REM OutputFile: "%2\%1.info"
REM ErrorFile: "%2\%1.err"

DEL "%2\%1.post.res"

REM delete the line before and uncomment the following line 
REM to execute the program
MOVE "%2\%1.dat" "%2\%1.m"

