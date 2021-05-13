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

    add_custom_command(
        OUTPUT
            ${HTML_DOCUMENTATIONS}/index.html
        DEPENDS
            neci_dev.md neci_user.md # TODO(@Oskar) Add fortran sources as dependencies
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${HTML_DOCUMENTATIONS} # Create output directory
        COMMAND
            ${FORD} ${PROJECT_SOURCE_DIR}/gen_docs.md --output_dir ${HTML_DOCUMENTATIONS}
        WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}
        COMMENT
            "Generating FORD documentation."
        VERBATIM
    )

    add_custom_target(doc_html
        DEPENDS
            ${HTML_DOCUMENTATIONS}/index.html
    )
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
            neci_dev.md neci_user.md packages_to_include.tex
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${PDF_DOCUMENTATIONS} # Create output directory
        COMMAND
            ${PROJECT_SOURCE_DIR}/tools/create_pdf_doc.sh neci_user.md packages_to_include.tex ${PDF_DOCUMENTATIONS}/neci_user.pdf
        COMMAND
            ${PROJECT_SOURCE_DIR}/tools/create_pdf_doc.sh neci_dev.md packages_to_include.tex ${PDF_DOCUMENTATIONS}/neci_dev.pdf
        WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}/docs
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
