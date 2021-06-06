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
    set(HTML_DOCUMENTATIONS ${PROJECT_BINARY_DIR}/documentation/html/)
    find_package(Python3 REQUIRED)
    find_program (FORD ford ${PROJECT_SOURCE_DIR}/External/ford/ford.py)

    get_absolute_sources(neci neci_sources)
    get_absolute_sources(libneci libneci_sources)
    set(sources ${neci_sources} ${libneci_sources})

    get_target_property(source_dir libneci SOURCE_DIR)

    add_custom_command(
        OUTPUT
            ${HTML_DOCUMENTATIONS}/index.html
        DEPENDS
            gen_docs.md pages/01_neci_user.md pages/02_neci_dev.md
            pages/tutorials/01_installation.md  pages/tutorials/02_basic_use.md
            pages/tutorials/FCIDUMP  pages/tutorials/index.md  pages/tutorials/n2.mol  pages/tutorials/plot_fcimcstats.plt
            "${sources}"
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${HTML_DOCUMENTATIONS} # Create output directory
        COMMAND
            ${Python3_EXECUTABLE} ${FORD} gen_docs.md --output_dir ${HTML_DOCUMENTATIONS} --src_dir=${source_dir}/lib
#             ${Python3_EXECUTABLE} ${FORD} gen_docs.md --output_dir ${HTML_DOCUMENTATIONS} --src_dir=${source_dir} --src_dir=${source_dir}/lib
        WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}/docs
        COMMENT
            "Generating FORD documentation."
        VERBATIM
    )

    add_custom_target(doc_html
        DEPENDS
            ${HTML_DOCUMENTATIONS}/index.html
    )
    # We rely on fypp preprocessed files by the normal build process
    add_dependencies(doc_html neci libneci)
    add_dependencies(doc doc_html)

endif()


find_program(PANDOC pandoc)
if (PANDOC)
    set(PDF_DOCUMENTATIONS ${PROJECT_BINARY_DIR}/documentation/pdf/)
    add_custom_command(
        OUTPUT
            ${PDF_DOCUMENTATIONS}/neci_user.pdf
            ${PDF_DOCUMENTATIONS}/neci_dev.pdf
        DEPENDS
            pages/01_neci_user.md pages/02_neci_dev.md pages/packages_to_include.tex ${PROJECT_SOURCE_DIR}/tools/create_pdf_doc.sh
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
