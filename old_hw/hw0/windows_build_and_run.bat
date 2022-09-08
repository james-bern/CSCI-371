@echo off
del executable.pdb
del executable.exe
cls

set OPTARG=d
IF NOT "%1"=="" set OPTARG=2
cl -O%OPTARG% ^
-W4 -wd4201 -wd4127 ^
/nologo -fp:except -Zi -GR- -EHa- -FC ^
/I.\ext\glfw\include ^
/I.\ext\glad\include ^
/EHsc /MDd hw.cpp ^
/Feexecutable.exe ^
/link /NODEFAULTLIB:MSVCRT ^
OpenGL32.lib user32.lib gdi32.lib shell32.lib vcruntime.lib ^
ext\glfw\windows\glfw3.lib

del executable.exp executable.ilk executable.lib hw.obj

start /i executable.exe
