add_executable(diff diff.c)
file(COPY_FILE "${CMAKE_CURRENT_SOURCE_DIR}/diff.dat" "${CMAKE_CURRENT_BINARY_DIR}/diff.dat")
