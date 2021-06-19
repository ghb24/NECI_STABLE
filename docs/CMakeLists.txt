# Add additional dependencies here.
# The order of files in each list matters!
# (It defines in which order they appear in the PDF)
set(markdown_files_user_doc
        pages/01_user_doc/index.md
        pages/01_user_doc/01_getting_into_the_game.md
        pages/01_user_doc/02_calculation_inputs.md
        pages/01_user_doc/03_calculation_outputs.md
        pages/01_user_doc/04_trial_wave_functions.md
        pages/01_user_doc/05_sampling_excited_states.md
        pages/01_user_doc/06_davidson_RAS.md
        pages/01_user_doc/07_error_analysis.md
        pages/01_user_doc/08_references.md
)

set(markdown_files_tutorial
        pages/02_tutorials/index.md
        pages/02_tutorials/01_installation.md
        pages/02_tutorials/02_basic_use.md
)

set(markdown_files_dev_doc
        pages/03_dev_doc/index.md
        pages/03_dev_doc/01_code_conventions.md
        pages/03_dev_doc/02_comments_and_documentation.md
        pages/03_dev_doc/03_review_guidelines.md
        pages/03_dev_doc/04_DRY.md
        pages/03_dev_doc/05_compiler_help.md
        pages/03_dev_doc/06_premature_optimization.md
        pages/03_dev_doc/07_memory_tracking.md
        pages/03_dev_doc/08_templating.md
        pages/03_dev_doc/09_testing.md
        pages/03_dev_doc/10_interfacing_C.md
        pages/03_dev_doc/11_debugging.md
        pages/03_dev_doc/12_important_globals.md
)

set(MEDIA_FILES
        media/N2_neci_files/FCIDUMP
        media/N2_neci_files/n2.mol
        media/N2_neci_files/plot_fcimcstats.plt
        media/plots1/blocking.png
        media/plots1/check_denominator.png
        media/plots1/check_numerator.png
        media/plots1/check_refPop.png
        media/plots1/check_shift_energy.png
        media/plots1/check_totE.png
        media/plots1/check_totWalkers.png
        media/plots2/blocking.png
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
    #
    # A/file1.txt A/file2.txt B/file1.txt -> A/ B/
    foreach(path ${${paths}})
        get_filename_component(dir ${path} DIRECTORY)
        if (NOT "${dir}" STREQUAL "")
            if (NOT ${dir} IN_LIST ${result})
                list(APPEND ${result} ${dir})
            endif()
        endif()
    endforeach()
    set(${result} ${${result}} PARENT_SCOPE)
endfunction()


find_package(Python3 REQUIRED)

set(PYTHON_VENV ${PROJECT_BINARY_DIR}/venv)
add_custom_command(
    OUTPUT
        ${PYTHON_VENV}
    COMMAND
        ${Python3_EXECUTABLE} -m venv ${PYTHON_VENV}
    COMMAND
        ${PYTHON_VENV}/bin/pip3 install ${PROJECT_SOURCE_DIR}/External/ford/
    COMMAND
        ${PYTHON_VENV}/bin/pip3 install bibtexparser
    COMMENT
        "Setting up virtual environment for ford."
)
set(custom_Python3 ${PYTHON_VENV}/bin/python3)


function(add_ford_documentation target configuration)
    # This function creates a target for the Ford documentation.
    # configuration can be one of neci;kneci;dneci;mneci;kdneci;kmneci
    #
    # The target is not built by default! Use add_dependencies if appropiate.

    if (NOT ${configuration} IN_LIST ${PROJECT_NAME}_CONFIGURATIONS)
        message(FATAL_ERROR "${configuration} has to be contained in [${${PROJECT_NAME}_CONFIGURATIONS}]")
    endif()

    set(THIS_HTML_DOC ${PROJECT_BINARY_DIR}/documentation/html/${configuration})

    find_program(FORD ford ${PROJECT_SOURCE_DIR}/External/ford/ford.py REQUIRED)
    # We want to install the custom ford into a virtual environment for the build process

    set(FORD_PAGES ${PROJECT_BINARY_DIR}/documentation/ford_pages/${configuration}/)

    set(markdown_files ${markdown_files_user_doc} ${markdown_files_tutorial} ${markdown_files_dev_doc})
    get_unique_dirs(markdown_files markdown_dirs)

    foreach(file ${markdown_files})
        get_filename_component(dir ${file} DIRECTORY)
        get_filename_component(dir ${dir} NAME)
        get_filename_component(name ${file} NAME)
        list(APPEND processed_files "${FORD_PAGES}/${dir}/${name}")
    endforeach()


    add_custom_command(
        OUTPUT
            ${processed_files}
        DEPENDS
            ${PYTHON_VENV}
            ${PROJECT_SOURCE_DIR}/tools/doc_compilation/process_markdown_for_ford.py
            pages/index.md
            ${markdown_files}
            LaTeX_headers/literature.bib
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${FORD_PAGES}
        COMMAND
            cp pages/index.md ${FORD_PAGES}
        COMMAND
            ${custom_Python3} ${PROJECT_SOURCE_DIR}/tools/doc_compilation/process_markdown_for_ford.py
                        --bibfile LaTeX_headers/literature.bib
                        --out_dir ${FORD_PAGES}
                        --md_dirs ${markdown_dirs}
                        --build_dir ${THIS_HTML_DOC}
        WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}/docs
        COMMENT
            "Preprocessing for Ford."
    )


    get_absolute_sources(${configuration} prog_sources)
    get_absolute_sources(lib${configuration} lib_sources)
    set(fortran_sources ${prog_sources} ${lib_sources})
    get_unique_dirs(fortran_sources source_dirs)
# TODO(Oskar, Philip): Put back all source files
#     foreach(dir ${source_dirs})
#         list(APPEND ford_options "--src_dir=${dir}")
#     endforeach()
    list(APPEND ford_options "--src_dir=${PROJECT_SOURCE_DIR}/src/lib")

    get_target_property(comp_definitions lib${configuration} COMPILE_DEFINITIONS)
    foreach(def ${comp_definitions})
        list(APPEND ford_options "--macro=${def}")
    endforeach()

    add_custom_command(
        OUTPUT
            ${THIS_HTML_DOC}/index.html
        DEPENDS
            gen_docs.md
            ${PYTHON_VENV}
            ${MEDIA_FILES}
            ${fortran_sources}
            ${processed_files}
            ${markdown_files}
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${THIS_HTML_DOC} # Create output directory
        COMMAND
            ${custom_Python3} ${FORD} gen_docs.md --output_dir ${THIS_HTML_DOC} ${ford_options} --page_dir ${FORD_PAGES}
        WORKING_DIRECTORY
            ${PROJECT_SOURCE_DIR}/docs
        COMMENT
            "Generating FORD documentation."
        VERBATIM
    )

    add_custom_target(${target}
        DEPENDS
            ${THIS_HTML_DOC}/index.html
    )
    # We rely on preprocessed (fypp and custom templating) files by the normal build process.
    # Ford will run the "normal" C-preprocessor itself.
    add_dependencies(${target} lib${configuration}_fypp_preprocessed lib${configuration}_templated)
endfunction()



add_custom_target(doc)

# We want our custom ford because of tailored changes!
# Do not look for a preinstalled ford, or at least make sure it works!
find_package(Git QUIET)
execute_process(
    COMMAND
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

    # Only build the documentation of `neci` by default (when typing `make doc`)
    add_dependencies(doc doc_html_neci)
endif()


find_program(PANDOC pandoc)
find_package(LATEX COMPONENTS PDFLATEX BIBER REQUIRED)
find_program(LATEXMK latexmk)
if (PANDOC AND LATEX_FOUND AND LATEXMK)
    # The PDF documentation is independent of possible configurations (neci, dneci etc.).
    set(PDF_DOCUMENTATIONS ${PROJECT_BINARY_DIR}/documentation/pdf/)
    set(TEX_FILES ${PROJECT_BINARY_DIR}/documentation/tex/)
    set(pandoc_options
            --include-in-header LaTeX_headers/packages_to_include.tex
            --standalone
            -V geometry:margin=3cm
            -V geometry:a4paper
            -V fontsize=18pt
            -V numbersections=true
            --listings
            --biblatex
                -V biblatexoptions:doi=true
                -V biblatexoptions:isbn=false
                -V biblatexoptions:url=false
                -V biblio-title=Bibliography
                --bibliography ${PROJECT_SOURCE_DIR}/docs/LaTeX_headers/literature.bib
            --toc
            -f markdown-tex_math_dollars+tex_math_single_backslash
    )


    function(add_tex_command markdown_files out_name comment)
        get_filename_component(out_dir ${out_name} DIRECTORY)
        add_custom_command(
            OUTPUT
                ${out_name}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/tools/doc_compilation/process_markdown_for_pdf.py
                ${PROJECT_SOURCE_DIR}/docs/LaTeX_headers/packages_to_include.tex
                ${${markdown_files}}
            COMMAND
                ${CMAKE_COMMAND} -E make_directory ${out_dir} # Create output directory
            COMMAND
                # the output of cat depends on the order in the lists!
                cat ${${markdown_files}}
                    | ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/doc_compilation/process_markdown_for_pdf.py
                    | ${PANDOC} ${pandoc_options} -o ${out_name} -
            WORKING_DIRECTORY
                ${PROJECT_SOURCE_DIR}/docs/
            COMMENT
                ${comment}
        )
    endfunction()

    function(add_pdf_command tex_file pdf_file comment)
        get_filename_component(tex_dir ${tex_file} DIRECTORY)
        get_filename_component(tex_stem ${tex_file} NAME_WE)
        get_filename_component(out_dir ${pdf_file} DIRECTORY)
        add_custom_command(
            OUTPUT
                ${pdf_file}
            DEPENDS
                ${tex_file}
                ${PROJECT_SOURCE_DIR}/docs/LaTeX_headers/packages_to_include.tex
                ${PROJECT_SOURCE_DIR}/docs/LaTeX_headers/literature.bib
            COMMAND
                ${CMAKE_COMMAND} -E make_directory ${out_dir} # Create output directory
            COMMAND
                ${LATEXMK} -view=pdf ${tex_file}
            COMMAND
                cp ${tex_stem}.pdf ${pdf_file}
            WORKING_DIRECTORY
                ${tex_dir}
            COMMENT
                ${comment}
        )
    endfunction()


    add_tex_command(markdown_files_user_doc ${TEX_FILES}/neci_user.tex "Generating user documentation tex file." )
    add_pdf_command(${TEX_FILES}/neci_user.tex ${PDF_DOCUMENTATIONS}/neci_user.pdf "Generating user pdf documentation.")
    add_tex_command(markdown_files_dev_doc ${TEX_FILES}/neci_dev.tex "Generating developer documentation tex file." )
    add_pdf_command(${TEX_FILES}/neci_dev.tex ${PDF_DOCUMENTATIONS}/neci_dev.pdf "Generating developer pdf documentation.")

    add_custom_target(doc_pdf
        DEPENDS
            ${PDF_DOCUMENTATIONS}/neci_user.pdf
            ${PDF_DOCUMENTATIONS}/neci_dev.pdf
    )
    add_dependencies(doc doc_pdf)
endif()
