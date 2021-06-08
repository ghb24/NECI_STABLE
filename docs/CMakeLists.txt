# Add additional dependencies here
set(MEDIA_FILES
        pages/tutorials/FCIDUMP
        pages/tutorials/n2.mol
        pages/tutorials/plot_fcimcstats.plt
        media/plots1/blocking.png media/plots1/check_denominator.png
        media/plots1/check_numerator.png
        media/plots1/check_refPop.png
        media/plots1/check_shift_energy.png
        media/plots1/check_totE.png
        media/plots1/check_totWalkers.png
        media/plots2/blocking.png
)


set(MANUALS_MARKDOWN
        pages/01_neci_user.md
        pages/02_neci_dev.md
)

set(ADDITIONAL_MARKDOWN_CONTENT
        pages/tutorials/index.md
        pages/tutorials/01_installation.md
        pages/tutorials/02_basic_use.md
)


function(get_absolute_sources target result)
    # This function returns the source files of a given target as absolute paths.
    get_target_property(rel_source_files ${target} SOURCES)
    get_target_property(source_dir ${target} SOURCE_DIR)

    foreach(path ${rel_source_files})
        if(NOT IS_ABSOLUTE ${path})
            set(path "${source_dir}/${path}")
        endif()
        list(APPEND ${result} ${path})
    endforeach()
    set(${result} ${${result}} PARENT_SCOPE)
endfunction()


function(get_unique_dirs paths result)
    # This function returns a unique list of all directories in paths.
    # paths have to be absolute.
    foreach(path ${${paths}})
        if(IS_ABSOLUTE ${path})
            get_filename_component(dir ${path} DIRECTORY)
            if (NOT ${dir} IN_LIST ${result})
                list(APPEND ${result} ${dir})
            endif()
        else()
            message(FATAL_ERROR "Non absolute path passed to `get_unique_dirs`")
        endif()
    endforeach()
    set(${result} ${${result}} PARENT_SCOPE)
endfunction()

function(add_ford_documentation target configuration)
    # This function creates a target for the Ford documentation.
    # configuration can be one of neci;kneci;dneci;mneci;kdneci;kmneci
    #
    # The target is not built by default! Use add_dependencies if appropiate.
    if (NOT ${configuration} IN_LIST ${PROJECT_NAME}_CONFIGURATIONS)
        message(FATAL_ERROR "${configuration} has to be contained in [${${PROJECT_NAME}_CONFIGURATIONS}]")
    endif()

    set(HTML_DOCUMENTATIONS ${PROJECT_BINARY_DIR}/documentation/html/${configuration})

    find_package(Python3 REQUIRED)
    find_program (FORD ford ${PROJECT_SOURCE_DIR}/External/ford/ford.py REQUIRED)

    get_absolute_sources(${configuration} prog_sources)
    get_absolute_sources(lib${configuration} lib_sources)
    set(fortran_sources ${prog_sources} ${lib_sources})
    get_unique_dirs(fortran_sources source_dirs)
    foreach(dir ${source_dirs})
        list(APPEND ford_options "--src_dir=${dir}")
    endforeach()

    get_target_property(comp_definitions lib${configuration} COMPILE_DEFINITIONS)
    foreach(def ${comp_definitions})
        list(APPEND ford_options "--macro=${def}")
    endforeach()

    add_custom_command(
        OUTPUT
            ${HTML_DOCUMENTATIONS}/index.html
        DEPENDS
            gen_docs.md
            ${MANUALS_MARKDOWN}
            ${ADDITIONAL_MARKDOWN_CONTENT}
            ${MEDIA_FILES}
            ${fortran_sources}
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${HTML_DOCUMENTATIONS} # Create output directory
        COMMAND
            ${Python3_EXECUTABLE} ${FORD} gen_docs.md --output_dir ${HTML_DOCUMENTATIONS} ${ford_options}
        WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}/docs
        COMMENT
            "Generating FORD documentation."
        VERBATIM
    )

    add_custom_target(${target}
        DEPENDS
            ${HTML_DOCUMENTATIONS}/index.html
    )
    # We rely on preprocessed files by the normal build process
    add_dependencies(${target} lib${configuration}_fypp_preprocessed lib${configuration}_templated)
endfunction()



add_custom_target(doc)

# We want our custom ford because of tailored changes!
# Do not look for a preinstalled ford, or at least make sure it works!
find_package(Git QUIET)
execute_process(COMMAND
                    ${GIT_EXECUTABLE} submodule update --init --recursive
                WORKING_DIRECTORY
                    ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE
                    GIT_SUBMOD_RESULT
)

if (GIT_SUBMOD_RESULT EQUAL "0")
    foreach(configuration ${${PROJECT_NAME}_CONFIGURATIONS})
        add_ford_documentation(doc_html_${configuration} ${configuration})
    endforeach()

    # Only build the documentation of `neci` by default
    add_dependencies(doc doc_html_neci)
endif()


find_program(PANDOC pandoc)
if (PANDOC)
    # The PDF documentation is independent of possible configurations.
    set(PDF_DOCUMENTATIONS ${PROJECT_BINARY_DIR}/documentation/pdf/)
    add_custom_command(
        OUTPUT
            ${PDF_DOCUMENTATIONS}/neci_user.pdf
            ${PDF_DOCUMENTATIONS}/neci_dev.pdf
        DEPENDS
            pages/packages_to_include.tex ${PROJECT_SOURCE_DIR}/tools/create_pdf_doc.sh
            ${MANUALS_MARKDOWN}
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${PDF_DOCUMENTATIONS} # Create output directory
        COMMAND
            ${PROJECT_SOURCE_DIR}/tools/create_pdf_doc.sh 01_neci_user.md packages_to_include.tex ${PDF_DOCUMENTATIONS}/neci_user.pdf
        COMMAND
            ${PROJECT_SOURCE_DIR}/tools/create_pdf_doc.sh 02_neci_dev.md packages_to_include.tex ${PDF_DOCUMENTATIONS}/neci_dev.pdf
        WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}/docs/pages/
        COMMENT
            "Generating user and developer PDF documentation."
        VERBATIM
    )

    add_custom_target(doc_pdf
        DEPENDS
            ${PDF_DOCUMENTATIONS}/neci_user.pdf
            ${PDF_DOCUMENTATIONS}/neci_dev.pdf
    )
    add_dependencies(doc doc_pdf)
endif()
