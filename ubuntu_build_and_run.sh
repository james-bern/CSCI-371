#!/bin/bash

if [ "$1" = "--help" ]; then
    echo "build and   run in   debug mode: [36m./ubuntu_build_and_run.sh[0m"
    echo "build and   run in release mode: [36m./ubuntu_build_and_run.sh --release[0m"
    echo "build and debug in      VS Code: [36m./ubuntu_build_and_run.sh --debug-vscode[0m"
else
    if [ -f "executable.exe" ]; then
        rm executable.exe
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

    g++ -Werror=vla -fno-strict-aliasing -O$OPTARG -g -std=c++11 -Wall -Wextra -Werror=vla -Wno-format-security -Wno-misleading-indentation -Wno-char-subscripts -Wno-conversion -Wno-write-strings -Wno-missing-field-initializers -o executable.exe main.cpp -lglfw -lGL

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
