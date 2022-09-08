rm -r executable.exe.dSYM
rm executable.exe

OPTARG=0
CONDITION=${1:-0}
if [ $# != 0 ]
then
OPTARG=3
fi

clang++ -O$OPTARG -std=c++11 -g -Wall -Wextra -Wno-missing-braces -Wno-missing-field-initializers -Wno-char-subscripts -Wno-writable-strings -I./ext/glfw/include -I./ext/glad/include -L./ext/glfw/mac -o executable.exe hw.cpp -lglfw3 -framework Cocoa -framework OpenGL -framework IOKit
./executable.exe
