@echo off

IF "%1"=="--help" (
    echo build and   run in   debug mode: [36mwindows_build_and_run.bat[0m
    echo build and   run in release mode: [36mwindows_build_and_run.bat --release[0m
    echo build and debug in      VS Code: [36mwindows_build_and_run.bat --debug-vscode[0m
    echo build and debug in     remedyBG: [36mwindows_build_and_run.bat --debug-remedybg[0m
) ELSE (
    IF EXIST "hw.obj"         ( del hw.obj         )
    IF EXIST "vc140.pdb"      ( del vc140.pdb      )
    IF EXIST "executable.pdb" ( del executable.pdb )
    IF EXIST "executable.exe" ( del executable.exe )
    IF EXIST "executable.exp" ( del executable.exp )
    IF EXIST "executable.ilk" ( del executable.ilk )
    IF EXIST "executable.lib" ( del executable.lib )

    cls

    (
        setlocal ENABLEDELAYEDEXPANSION
        set OPTARG=2
        set DEBARG=-Zi
        IF "%1"=="--release" (
            echo [36m[cow] compiling in release mode[0m
        ) ELSE (
            echo [36m[cow] compiling in debug mode[0m
            set OPTARG=d
        )

        (
            cl -O!OPTARG! ^
            -W4 -wd4201 -wd4127 ^
            /nologo -fp:except !DEBARG! -GR- -EHa- -FC ^
            /I.\ext\glfw\include ^
            /I.\ext\glad\include ^
            /EHsc /MDd hw.cpp ^
            /Feexecutable.exe ^
            /link /NODEFAULTLIB:MSVCRT ^
            OpenGL32.lib user32.lib gdi32.lib shell32.lib vcruntime.lib ^
            ext\glfw\windows\glfw3.lib
        )
    )

    call :setESC

    IF EXIST "executable.exe" (
        IF "%1"=="--debug-vscode" (
            echo [36m[cow] debugging in Visual Studio Code[0m
            _xplat_debug_vscode.bat
        ) ELSE IF "%1"=="--debug-remedybg" (
            echo [36m[cow] debugging in remedyBG[0m
            call _windows_debug_remedybg.bat
        ) ELSE (
            echo [36m[cow] running executable[0m
            @echo on
            start executable.exe
        )
    )
)
@echo off
:setESC
for /F %%a in ('"prompt $E$S & echo on & for %%b in (1) do rem"') do set "ESC=%%a"
exit /B


