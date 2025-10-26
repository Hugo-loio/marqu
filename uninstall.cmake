set(MANIFEST "${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt")

if(NOT EXISTS "${MANIFEST}")
    message(FATAL_ERROR "Cannot find install manifest: ${MANIFEST}")
endif()

file(READ "${MANIFEST}" files)

# Split the contents into a list of files, removing extra spaces or newlines
string(REPLACE "\n" ";" files "${files}")

foreach(file ${files})
    message(STATUS "Removing: ${file}")
    file(REMOVE ${file})
endforeach()
