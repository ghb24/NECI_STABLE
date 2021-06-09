# Add additional dependencies here.
# The order of files in each list matters!
# (It defines in which order they appear in the PDF)
set(markdown_files_user_doc
        pages/01_user_doc/index.md
        pages/01_user_doc/01_getting_into_the_game.md
        pages/01_user_doc/02_calculation_inputs.md
        pages/01_user_doc/03_trial_wave_functions.md
        pages/01_user_doc/04_sampling_excited_states.md
        pages/01_user_doc/05_davidson_RAS.md
        pages/01_user_doc/06_error_analysis.md
        pages/01_user_doc/07_references.md
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

set(markdown_files_tutorial
    pages/02_tutorials/index.md
    pages/02_tutorials/01_installation.md
    pages/02_tutorials/02_basic_use.md
)

set(MEDIA_FILES
        pages/02_tutorials/FCIDUMP
        pages/02_tutorials/n2.mol
        pages/02_tutorials/plot_fcimcstats.plt
        media/plots1/blocking.png media/plots1/check_denominator.png
        media/plots1/check_numerator.png
        media/plots1/check_refPop.png
        media/plots1/check_shift_energy.png
        media/plots1/check_totE.png
        media/plots1/check_totWalkers.png
        media/plots2/blocking.png
)

find_package(Python3 REQUIRED)

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
    set(PYTHON_VENV ${PROJECT_BINARY_DIR}/venv)

    find_program(FORD ford ${PROJECT_SOURCE_DIR}/External/ford/ford.py REQUIRED)
    # We want to install the custom ford into a virtual environment for the build process
    add_custom_command(
        OUTPUT
            ${PYTHON_VENV}
        COMMAND
            ${Python3_EXECUTABLE} -m venv ${PYTHON_VENV}
        COMMAND
            ${PYTHON_VENV}/bin/pip3 install ${PROJECT_SOURCE_DIR}/External/ford/
        COMMENT
            "Setting up virtual environment for ford."
    )
    set(custom_Python3 ${PYTHON_VENV}/bin/python3)


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
            ${PYTHON_VENV}
            ${markdown_files_user_doc}
            ${markdown_files_tutorial}
            ${markdown_files_dev_doc}
            ${MEDIA_FILES}
            ${fortran_sources}
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${HTML_DOCUMENTATIONS} # Create output directory
        COMMAND
            ${custom_Python3} ${FORD} gen_docs.md --output_dir ${HTML_DOCUMENTATIONS} ${ford_options}
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
if (PANDOC)
    # The PDF documentation is independent of possible configurations (neci, dneci etc.).
    set(PDF_DOCUMENTATIONS ${PROJECT_BINARY_DIR}/documentation/pdf/)
    set(pandoc_options
            --listings
            --include-in-header LaTeX_headers/packages_to_include.tex
            -V geometry:margin=3cm
            -V geometry:a4paper
            -V fontsize=18pt
            -V numbersections=true
            --toc
            -f markdown-tex_math_dollars+tex_math_single_backslash
    )


    function(add_pdf_command markdown_files out_name comment)
        add_custom_command(
            OUTPUT
                ${PDF_DOCUMENTATIONS}/${out_name}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/tools/create_pdf_doc.py
                LaTeX_headers/packages_to_include.tex
                ${${markdown_files}}
            COMMAND
                ${CMAKE_COMMAND} -E make_directory ${PDF_DOCUMENTATIONS} # Create output directory
            COMMAND
                cat ${${markdown_files}}
                | ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/create_pdf_doc.py
                | tee ~/Desktop/test.md
                | ${PANDOC} ${pandoc_options} -o ${PDF_DOCUMENTATIONS}/${out_name} -
            WORKING_DIRECTORY
                ${PROJECT_SOURCE_DIR}/docs/
            COMMENT
                ${comment}
        )
    endfunction()

    add_pdf_command(markdown_files_user_doc neci_user.pdf "Generating user PDF documentation." )
    add_pdf_command(markdown_files_dev_doc neci_dev.pdf "Generating developer PDF documentation." )

    add_custom_target(doc_pdf
        DEPENDS
            ${PDF_DOCUMENTATIONS}/neci_user.pdf
            ${PDF_DOCUMENTATIONS}/neci_dev.pdf
    )
    add_dependencies(doc doc_pdf)
endif()
