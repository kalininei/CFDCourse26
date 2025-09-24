#!/bin/bash
CONTAINER_NAME="cfd26"

# Функция для поиска корня проекта
find_project_root() {
    local dir=$(pwd)
    while [ "$dir" != "/" ]; do
        # Ищем маркеры корня проекта
        if [ -d "$dir/.git" ] || [ -f "$dir/.projectroot" ]; then
            echo "$dir"
            return 0
        fi
        dir=$(dirname "$dir")
    done
    # Если корень не найден, используем текущую директорию
    echo $(pwd)
}

# Находим корень проекта на хосте
HOST_PROJECT_ROOT=$(find_project_root)

# Определяем относительный путь от корня до текущей директории
RELATIVE_PATH=$(realpath --relative-to="$HOST_PROJECT_ROOT" $(pwd))

# В контейнере: /app + относительный путь
CONTAINER_WORK_DIR="/app/$RELATIVE_PATH"

echo "Host project root: $HOST_PROJECT_ROOT" >&2
echo "Relative path: $RELATIVE_PATH" >&2
echo "Container work dir: $CONTAINER_WORK_DIR" >&2

# Запускаем компилятор в контейнере
docker exec -i -w "$CONTAINER_WORK_DIR" "$CONTAINER_NAME" gcc "$@"
