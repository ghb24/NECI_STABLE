.. _theory_introduction:

------------
Introduction
------------

The energy of a system can be evaluated using a standard statisical mechanics result:

.. math::
 E = \frac{\operatorname{Tr}[H e^{-\beta H}]}{\operatorname{Tr}[e^{-\beta H}]}

We choose to work in a Slater Determinant space, which is, by construction, anti-symmetric.  In this space the energy expression becomes:

.. math::
 E &= \frac{\sum_{\veci} \bra D_{\veci} | H e^{-\beta H} | D_{\veci} \ket}{\sum_{veci} \bra D_{\veci} | e^{-\beta H} | D_{\veci} \ket} \\

   &= \frac{\sum_{\veci} w_{\veci} \tilde{E}_{\veci}}{\sum_{\veci} w_{\veci}}

A given term in the numerator is simply the differential of the
corresponding term in the denominator.  There is a cleaner and more
efficient way of evaluating the numerator than differentiation, but we
will first turn our attention to the denominator.

We can expand each term into a closed path of :math:`P` steps through the discrete Slater Determinant space:

.. math::
     w^P_{\veci,\veci_1,\veci_2,\cdots,\veci_p,\veci} &= \sum_{\veci_1} \sum_{\veci_1} \cdots \sum_{\veci_P} \bra D_{\veci_1} | e^{\beta H/P} | D_{\veci_2}  \ket \bra D_{\veci_2} | e^{\beta H/P} | D_{\veci_3}  \ket \cdots \bra D_{\veci_P} | e^{\beta H/P} | D_{\veci_1}  \ket \\

     & =  \sum_{\veci_1} \sum_{\veci_1} \cdots \sum_{\veci_P} \rho_{\veci_1\veci_2} \rho_{\veci_2\veci_3} \cdots \rho_{\veci_P\veci_1},

where the :math:`\rho` matrix consists of elements
:math:`\rho_{\veci\vecj}=\bra D_{\veci} | e^{\beta H/P} | D_{\vecj} \ket`.

Each path does not necessarily visit :math:`P-1` determinants: "hopping"
terms are allowed.  Due to the :math:`\rho` matrix being diagonally
dominant, paths containing small numbers of unique determinants will
tend to have a much greater contribution to the overall energy.

The size of the Slater determinant space grows factorially with the number
of electrons and virtual orbitals, making it impossible to sum together
all the paths.  Furthermore, the sign of a path is an incredibly poorly
behaved quantity.  It is possible to perform an analytical resummation
of the paths into objects we term graphs, where each graph contains
paths which only visit the vertices contained within the graph.

.. note::
    To come: pictures of paths --> graph.

The expression for the energy now becomes a sum over Slater determinants and a sum graphs which originate from each Slater determinant:

.. math::
  E = \frac{\sum_{\veci} \sum_G w_{\veci}[G] \tilde{E}_{\veci}[G]}{\sum_{\veci} \sum_G w_{\veci}[G]}

Furthermore, the graphs have a much better sign behaviour: there are many graphs with a definite-positive weight, at least for graphs with less than 5 vertices, which makes a Monte Carlo approach feasible.

The resummation of paths into graphs still leaves a sum that is far too
large to be completely evaluated.  There are various approximations we
can apply.

     #. Use a single reference reference approach, i.e. approximate the 
        energy with:

        .. math::
           E = \frac{\sum_G w_{\vecz}[G] \tilde{E}_{\vecz}[G]}{\sum_G w_{\vecz}[G]}

        where :math:`\vecz` refers to the reference (i.e. Hartree--Fock)
        determinant.

        This sum, in general, still contains too many terms (and grows
        too rapidly with system size) to be of much use.
        
     #. Truncate the sum at a certain graph size (e.g. restrict it to
        two or three vertices).  This approach is referred to as a
        **VERTEX SUM** approach.
     #. Find a large graph that is a good approximation to the ground state
        and is easy to evaluate.  Our current model is the single and
        double excitation star, which contains all single and double
        excitations connected to the reference determinant but ignores any
        connections between the excited determinants.  In other words,
        it couples all the single and double excitations together,
        but only through the reference determinant.  This method is
        referred to as a **VERTEX STAR** approach.  The star contains all
        graphs in the sum truncated at the two vertex level and much more but
        at no additional costly integrals to evaluate.  This makes it a very
        attractive approach.
        
.. note::
    To come: the propogation operator.
