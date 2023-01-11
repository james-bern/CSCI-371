@echo off
setlocal ENABLEDELAYEDEXPANSION

Taskkill /IM "executable.exe" /F  >nul 2>&1

for %%a in (%*) do set "argv[%%a]=1"

IF defined argv[--help] (
    echo build and   run in   debug mode: [36mwindows_build_and_run.bat[0m
    echo build and   run in release mode: [36mwindows_build_and_run.bat --release[0m
    echo build and debug in      VS Code: [36mwindows_build_and_run.bat --debug-vscode[0m
    echo build and debug in     remedyBG: [36mwindows_build_and_run.bat --debug-remedybg[0m
    echo ---
    echo include [36m--eigen[0m to use Eigen's sparse linear solver
) ELSE (
    IF EXIST "main.obj"       ( del hw.obj         )
    IF EXIST "vc140.pdb"      ( del vc140.pdb      )
    IF EXIST "executable.pdb" ( del executable.pdb )
    IF EXIST "executable.exe" ( del executable.exe )
    IF EXIST "executable.exp" ( del executable.exp )
    IF EXIST "executable.ilk" ( del executable.ilk )
    IF EXIST "executable.lib" ( del executable.lib )

    cls

    (
    set OPTARG=2
    set DEBARG=-Zi
    IF defined argv[--release] (
        echo [36m[cow] compiling in release mode[0m
    ) ELSE (
        echo [36m[cow] compiling in debug mode[0m
        set OPTARG=d
    )

    set EIGEN_LIB=
    set EIGEN_DEFINE=
    IF defined argv[--eigen] (
        echo [36m[cow] linking against precompiled Eigen linear solver[0m
        set EIGEN_LIB=codebase\ext\windows_eigen.lib
        set EIGEN_DEFINE=/DUSE_EIGEN
    )

    (
    cl -O!OPTARG! ^
    -W4 -wd4201 -wd4127 ^
    /nologo -fp:except !DEBARG! -GR- -EHa- -FC ^
    /I.\codebase\ext\ ^
    !EIGEN_DEFINE! ^
    /EHsc /MDd main.cpp ^
    /Feexecutable.exe ^
    /link /NODEFAULTLIB:MSVCRT ^
    OpenGL32.lib user32.lib gdi32.lib shell32.lib vcruntime.lib ^
    codebase\ext\windows_glfw3.lib ^
    !EIGEN_LIB!
)
    )

    call :setESC

    IF EXIST "executable.exe" (
        IF defined argv[--debug-vscode] (
            echo [36m[cow] debugging in Visual Studio Code[0m
            _xplat_debug_vscode.bat
        ) ELSE IF defined argv[--debug-remedybg] (
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
endlocal

IF EXIST "main.obj"       ( del main.obj       )
IF EXIST "executable.ilk" ( del executable.ilk )

:setESC
for /F %%a in ('"prompt $E$S & echo on & for %%b in (1) do rem"') do set "ESC=%%a"
exit /B


