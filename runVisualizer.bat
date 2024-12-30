
@echo off
rem Set the path to the Rscript executable
set RSCRIPT="C:/PROGRA~1/R/R-44~1.1/bin/x64/R.exe"
rem Set the path to the R script to execute
set RSCRIPT_FILE="C:\Users\user\Load_All_Models.R"
rem set RSCRIPT_FILE="C:\Users\user\Rice_Model.R"
rem set RSCRIPT_FILE="C:\Users\user\RW_SystemTool.rmd"
rem set RSCRIPT_FILE="C:\Users\user\Wheat_Model.R"
rem set RSCRIPT_FILE="C:\Users\user\Rice_Wheat_Model.R"
rem Execute the R script
%RSCRIPT% %RSCRIPT_FILE%
rem Pause so the user can see the output
pause