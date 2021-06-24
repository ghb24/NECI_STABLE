# This is the main target of this file
add_custom_target(doc)

# Add additional dependencies here.
# The order of files in each list matters!
# (It defines in which order they appear in the PDF)
set(markdown_files_user_doc
        pages/01_user_doc/index.md
        pages/01_user_doc/01_getting_into_the_game.md
        pages/01_user_doc/02_intro_tutorial.md
        pages/01_user_doc/03_calculation_inputs.md
        pages/01_user_doc/03_calculation_outputs.md
        pages/01_user_doc/04_trial_wave_functions.md
        pages/01_user_doc/05_sampling_excited_states.md
        pages/01_user_doc/06_davidson_RAS.md
        pages/01_user_doc/07_error_analysis.md
        pages/01_user_doc/08_pylib.md
        pages/01_user_doc/09_references.md
)

set(markdown_files_dev_doc
        pages/02_dev_doc/index.md
        pages/02_dev_doc/01_code_conventions.md
        pages/02_dev_doc/02_comments_and_documentation.md
        pages/02_dev_doc/03_git_overview.md
        pages/02_dev_doc/04_review_guidelines.md
        pages/02_dev_doc/05_DRY.md
        pages/02_dev_doc/06_compiler_help.md
        pages/02_dev_doc/06_premature_optimization.md
        pages/02_dev_doc/07_memory_tracking.md
        pages/02_dev_doc/08_templating.md
        pages/02_dev_doc/09_testing.md
        pages/02_dev_doc/10_interfacing_C.md
        pages/02_dev_doc/11_debugging.md
        pages/02_dev_doc/12_important_globals.md
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

set(DOC_ROOT_SRC ${PROJECT_SOURCE_DIR}/docs)

set(DOC_WORKDIR ${PROJECT_BINARY_DIR}/documentation/work_dir)
set(MEDIA_WORKDIR ${DOC_WORKDIR}/media)
set(TEX_WORKDIR ${DOC_WORKDIR}/tex)
set(FORD_ROOT_WORKDIR ${DOC_WORKDIR}/ford)


foreach(file ${MEDIA_FILES})
    list(APPEND processed_media_files "${DOC_WORKDIR}/${file}")
endforeach()

add_custom_command(
    OUTPUT
        ${processed_media_files}
    DEPENDS
        ${MEDIA_FILES}
    COMMAND
        ${CMAKE_COMMAND} -E make_directory ${DOC_WORKDIR} # Create output directory
    COMMAND
        cp -r media ${MEDIA_WORKDIR}
    WORKING_DIRECTORY
        ${PROJECT_SOURCE_DIR}/docs
    COMMENT
        "Copying media files"
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

find_package(Git QUIET REQUIRED)
execute_process(
    COMMAND
        ${GIT_EXECUTABLE} submodule update --init --recursive
    WORKING_DIRECTORY
        ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE
        GIT_SUBMOD_RESULT
)

if (NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(FATAL_ERROR "Git submodules required. Check your internet connection.")
endif()


find_package(Python3 REQUIRED)
set(PYTHON_VENV ${PROJECT_BINARY_DIR}/venv)
set(custom_Python3 ${PYTHON_VENV}/bin/python3)
add_custom_command(
    OUTPUT
        ${PYTHON_VENV}
    COMMAND
        ${Python3_EXECUTABLE} -m venv ${PYTHON_VENV}
    COMMAND
        ${custom_Python3} -m pip install wheel --no-warn-script-location --quiet
    COMMENT
        "Setting up virtual python environment."
)


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

    set(FORD_WORKDIR ${FORD_ROOT_WORKDIR}/${configuration})
    set(FORD_PAGES ${FORD_WORKDIR}/pages)

    set(markdown_files ${markdown_files_user_doc} ${markdown_files_dev_doc})
    get_unique_dirs(markdown_files markdown_dirs)

    foreach(path ${markdown_files})
        list(APPEND src_files "${DOC_ROOT_SRC}/${path}")
        get_filename_component(dir ${path} DIRECTORY)
        get_filename_component(dir ${dir} NAME)
        get_filename_component(name ${path} NAME)
        list(APPEND processed_files "${FORD_PAGES}/${dir}/${name}")
    endforeach()


    add_custom_command(
        OUTPUT
            ${processed_files}
        DEPENDS
            ${PYTHON_VENV}
            ${PROJECT_SOURCE_DIR}/tools/doc_compilation/process_markdown_for_ford.py
            gen_docs.md
            pages/index.md
            ${src_files}
            LaTeX-headers/literature.bib
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${FORD_PAGES}
        COMMAND
            ${custom_Python3} -m pip install bibtexparser --quiet
        COMMAND
            cp gen_docs.md ${FORD_WORKDIR}
        COMMAND
            cp pages/index.md ${FORD_PAGES}
        COMMAND
            ${custom_Python3} ${PROJECT_SOURCE_DIR}/tools/doc_compilation/process_markdown_for_ford.py
                        --bibfile LaTeX-headers/literature.bib
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
    foreach(dir ${source_dirs})
        list(APPEND ford_options "--src_dir=${dir}")
    endforeach()
# Use the following instead if you want to compile ford fast
#     list(APPEND ford_options "--src_dir=${PROJECT_SOURCE_DIR}/src/lib")

    get_target_property(include_dirs lib${configuration} INCLUDE_DIRECTORIES)
    foreach(dir ${include_dirs})
        list(APPEND ford_options "--include=${dir}")
    endforeach()

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
            ${processed_media_files}
            ${fortran_sources}
            ${processed_files}
            ${markdown_files}
        # We want to install the custom ford into a virtual environment for the build process
        # We want our custom ford because of tailored changes!
        # Do not look for a preinstalled ford, or at least make sure it works!
        COMMAND
            ${custom_Python3} -m pip install ${PROJECT_SOURCE_DIR}/External/ford/ --quiet
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${THIS_HTML_DOC} # Create output directory
        COMMAND
            ${custom_Python3} ${FORD} gen_docs.md --output_dir ${THIS_HTML_DOC} ${ford_options} --page_dir ${FORD_PAGES}
        WORKING_DIRECTORY
            ${FORD_WORKDIR}
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



foreach(configuration ${${PROJECT_NAME}_CONFIGURATIONS})
    add_ford_documentation(doc_html_${configuration} ${configuration})
endforeach()

# Only build the documentation of `neci` by default (when typing `make doc`)
add_dependencies(doc doc_html_neci)


find_program(PANDOC pandoc REQUIRED)
find_package(LATEX COMPONENTS PDFLATEX BIBER REQUIRED)
find_program(LATEXMK latexmk REQUIRED)
if (PANDOC AND LATEX_FOUND AND LATEXMK)
    # The PDF documentation is independent of possible configurations (neci, dneci etc.).
    set(PDF_DOCUMENTATIONS ${PROJECT_BINARY_DIR}/documentation/pdf/)
    set(pandoc_options
            --include-in-header LaTeX-headers/packages_to_include.tex
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
                --bibliography literature.bib
            --toc
            -f markdown-tex_math_dollars+tex_math_single_backslash
    )


    function(add_tex_command markdown_files out_name comment)
        foreach(path ${${markdown_files}})
            list(APPEND src_files "${DOC_ROOT_SRC}/${path}")
        endforeach()

        get_filename_component(out_dir ${out_name} DIRECTORY)
        add_custom_command(
            OUTPUT
                ${out_name} ${out_dir}/literature.bib
            DEPENDS
                ${PROJECT_SOURCE_DIR}/tools/doc_compilation/process_markdown_for_pdf.py
                ${PROJECT_SOURCE_DIR}/docs/LaTeX-headers/packages_to_include.tex
                ${src_files}
            COMMAND
                ${CMAKE_COMMAND} -E make_directory ${out_dir} # Create output directory
            COMMAND
                # the output of cat depends on the order in the lists!
                cat ${${markdown_files}}
                    | ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/doc_compilation/process_markdown_for_pdf.py --media_dir ${MEDIA_WORKDIR}
                    | ${PANDOC} ${pandoc_options} -o ${out_name} -
            COMMAND
                cp LaTeX-headers/literature.bib ${out_dir}/literature.bib
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
                ${PROJECT_SOURCE_DIR}/docs/LaTeX-headers/packages_to_include.tex
                ${tex_dir}/literature.bib
                ${processed_media_files}
            COMMAND
                ${CMAKE_COMMAND} -E make_directory ${out_dir} # Create output directory
            COMMAND
                ${LATEXMK} -view=pdf ${tex_file} 1> /dev/null
            COMMAND
                cp ${tex_stem}.pdf ${pdf_file}
            WORKING_DIRECTORY
                ${tex_dir}
            COMMENT
                ${comment}
        )
    endfunction()


    add_tex_command(markdown_files_user_doc ${TEX_WORKDIR}/user/neci_user.tex "Generating user documentation tex file." )
    add_pdf_command(${TEX_WORKDIR}/user/neci_user.tex ${PDF_DOCUMENTATIONS}/neci_user.pdf "Generating user pdf documentation.")
    add_tex_command(markdown_files_dev_doc ${TEX_WORKDIR}/dev/neci_dev.tex "Generating developer documentation tex file." )
    add_pdf_command(${TEX_WORKDIR}/dev/neci_dev.tex ${PDF_DOCUMENTATIONS}/neci_dev.pdf "Generating developer pdf documentation.")

    add_custom_target(doc_pdf
        DEPENDS
            ${PDF_DOCUMENTATIONS}/neci_user.pdf
            ${PDF_DOCUMENTATIONS}/neci_dev.pdf
    )
    add_dependencies(doc doc_pdf)
endif()
