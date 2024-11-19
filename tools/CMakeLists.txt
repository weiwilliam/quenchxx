# Install scripts
list( APPEND tool_files
    check_diff.sh
    check_diff_recursively.sh
    cpplint.py
)

foreach(FILENAME IN LISTS tool_files)
    set(_tool ${PROJECT_NAME}_${FILENAME})
    #Install in build-tree
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
        ${CMAKE_BINARY_DIR}/bin/${_tool} )
    #Install in install-tree
    install(PROGRAMS ${FILENAME} RENAME ${_tool} DESTINATION bin )

    #TODO remove this once no packages depend on these tools in their original location.
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
        ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach()
