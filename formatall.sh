find src/ -type f \( -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.hpp" -o -name "*.cc" \) -print0 | xargs -0 clang-format -i -style=file
