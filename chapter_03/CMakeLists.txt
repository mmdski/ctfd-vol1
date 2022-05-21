add_executable(diff diff.c)
file(COPY_FILE "${CMAKE_CURRENT_SOURCE_DIR}/diff.dat" "${CMAKE_CURRENT_BINARY_DIR}/diff.dat")

add_executable(prob3.7 prob3.7.c)
file(COPY_FILE "${CMAKE_CURRENT_SOURCE_DIR}/prob3.7.dat" "${CMAKE_CURRENT_BINARY_DIR}/prob3.7.dat")
