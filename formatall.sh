#!/bin/bash
set -e

exclude_paths=(-not -path "./.*" -not -path "./build*" -not -path "*/venv/*" -not -path "*/env/*" -not -path "*/__pycache__/*")

echo "=== START FORMATTING ==="

format_file() {
    local formatter="$1"
    local file="$2"

    local original_mtime=$(stat -c %Y "$file")
    local original_content=$(cat "$file")

    $formatter "$file" 2>/dev/null || true

    local new_content=$(cat "$file")
    if [ "$original_content" = "$new_content" ]; then
        touch -d "@$original_mtime" "$file"
    else
	echo " format: $file"
    fi
}

echo "CHECK C++ ..."
find . \( -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -o -name "*.cc" -o -name "*.cxx" \) "${exclude_paths[@]}" | \
while read file; do
    format_file "clang-format -i -style=file" "$file"
done

echo "CHECK PYTHON..."
find . \( -name "*.py" \) "${exclude_paths[@]}" | \
while read file; do
    format_file "black -q" "$file"
done

echo "CHECK CMake..."
find . \( -name "CMakeLists.txt" -o -name "*.cmake" \) "${exclude_paths[@]}" | \
while read file; do
    format_file "cmake-format -i" "$file"
done

echo "=== DONE ==="
