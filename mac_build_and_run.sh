if [ "$1" = "--help" ]; then
    echo "build and   run in   debug mode: [36m./mac_build_and_run.bat[0m"
    echo "build and   run in release mode: [36m./mac_build_and_run.bat --release[0m"
    echo "build and debug in      VS Code: [36m./mac_build_and_run.bat --debug-vscode[0m"
else
    if [ -f "executable.exe" ]; then
        rm executable.exe
        rm -r executable.exe.dSYM
    fi

    clear

    OPTARG=0
    CONDITION=${1:-0}
    if [ "$1" = "--release" ]; then
        echo "[36m[cow] compiling in release mode[0m"
        OPTARG=3
    else
        echo "[36m[cow] compiling in debug mode[0m"
    fi

    clang++ -Wshadow -Werror=vla -fno-strict-aliasing -O$OPTARG -std=c++11 -g -Wall -Wextra -Wno-deprecated-declarations -Wno-missing-braces -Wno-missing-field-initializers -Wno-char-subscripts -Wno-writable-strings -I./codebase/ext -I./codebase/ext -L./codebase/ext -o executable.exe main.cpp -lglfw3 -framework Cocoa -framework OpenGL -framework IOKit -framework AudioUnit -framework CoreAudio -framework AudioToolbox

    if [ -f "executable.exe" ]; then
        if [ "$1" = "--debug-vscode" ]; then
            echo "[36m[cow] debugging in Visual Studio Code[0m"
            source _xplat_debug_vscode.bat
        else
            echo "[36m[cow] running executable[0m"
            ./executable.exe
        fi
    fi
fi
