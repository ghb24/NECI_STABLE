macro( neci_print_summary )


message(STATUS "---------------------------------------------------------")
message(STATUS "Target libraries: ${${PROJECT_NAME}_ALL_LIBS}")
message(STATUS "Target executables: ${${PROJECT_NAME}_ALL_EXES}")


message(STATUS "")
message(STATUS "  +---------------------------------------------------+")
message(STATUS "  | ${PROJECT_NAME} configuration now complete                   |")
message(STATUS "  |                                                   |")
message(STATUS "  | VERSION ${${PROJECT_NAME}_VERSION_STR}                                     |")
message(STATUS "  | SHAID   ${${PROJECT_NAME}_GIT_SHA1}  |")
message(STATUS "  +---------------------------------------------------+")
message(STATUS "")
message(STATUS "  You can now do 'make' to compile the software.")
message(STATUS "")

endmacro()
