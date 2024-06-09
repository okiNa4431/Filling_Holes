set(BUILD_TARGET "Filling_Holes")

#実行ファイル生成
add_executable(${BUILD_TARGET})
target_sources(${BUILD_TARGET} PUBLIC "main.cpp" "Filling_Holes.h" "Filling_Holes.cpp" ${COMMON_HEADERS})