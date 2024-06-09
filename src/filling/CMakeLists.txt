set(BUILD_TARGET "Filling_Holes")

#実行ファイル生成
add_executable(${BUILD_TARGET})
target_sources(${BUILD_TARGET} PUBLIC "main.cpp" "Filling_Holes.h" "Filling_Holes.cpp" ${COMMON_HEADERS})

#Eigenの設定
target_include_directories(${BUILD_TARGET} PUBLIC /Libraries/eigen-3.4.0)