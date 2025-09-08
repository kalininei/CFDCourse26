#!/bin/bash
set -e

exclude_paths=(-not -path "./.*" -not -path "./build*" -not -path "*/venv/*" -not -path "*/env/*" -not -path "*/__pycache__/*")

echo "=== START FORMATTING ==="

format_file() {
    local formatter="$1"
    local file="$2"
    local extension="${file##*.}"
    
    local temp_file=$(mktemp)
    cp "$file" "$temp_file"
    
    $formatter "$temp_file" >/dev/null 2>&1
    
    if ! diff -q "$file" "$temp_file" >/dev/null; then
        echo " format: $file"
        mv "$temp_file" "$file"
    else
        rm -f "$temp_file"
    fi
}

# 1. C++ файлы
echo "CHECK C++ ..."
find . \( -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -o -name "*.cc" -o -name "*.cxx" \) "${exclude_paths[@]}" | \
while read file; do
    format_file "clang-format -i -style=file" "$file"
done

# 2. Python файлы  
echo "CHECK PYTHON..."
find . \( -name "*.py" \) "${exclude_paths[@]}" | \
while read file; do
    format_file "black -q" "$file"
done

# 3. CMake файлы
echo "CHECK CMake..."
find . \( -name "CMakeLists.txt" -o -name "*.cmake" \) "${exclude_paths[@]}" | \
while read file; do
    format_file "cmake-format -i" "$file"
done

echo "=== DONE ==="
