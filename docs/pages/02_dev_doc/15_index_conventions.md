---
title: Indexing conventions
---

# Index Conventions

There are two main conventions in the quantum-chemistry community to
index two-electron integrals which also determines how to write the matrix elements.

For a given first quantised two-electron operator \(g^c(x_1, x_2)\),
we can define integrals over the orbitals \(\phi\)
using the so called **chemist's notation** \(g_{PQRS}\) or
the **physicist's notation** \(U_{PRQS}\) by:
\begin{equation}
    g_{PQRS}
    = U_{PRQS}
    =
        \int \int
        \phi^*_{P}(x_1)
        \phi^*_{R}(x_2)
        g^c(x_1, x_2)
        \phi_{Q}(x_1)
        \phi_{S}(x_2)
        \,
        \mathrm{d} x_1
        \mathrm{d} x_2
    \quad.
\end{equation}

With these integrals we can write the second-quantised
two electron operator in both notations as:
\begin{equation}
        \hat{g}
    =
        \sum_{PQRS}
        a^\dagger_P
        a^\dagger_R
        a_S
        a_Q
        g_{PQRS}
    =
        \sum_{PQRS}
        a^\dagger_P
        a^\dagger_R
        a_S
        a_Q
        U_{PRQS}
\end{equation}


Typical textbooks that assume the chemist's notation are the purple book[@PurpleBook] or Szabo-Ostlund[@Szabo].

The `FCIDUMP` file that transfers the electronic integrals from
`Molpro` or `Molcas` to `NECI` assumes the **chemist's notation**.

Internally `NECI` uses the **physicist's notation**,
i.e. the function `get_umat_el` that returns the stored two-electronic integrals
uses the indexing of \(U_{PRQS}\).
(Which is also why it is called `umat` in the first place.)

If we write a double excitation with second quantised operators,
we have the following commutation relation:
\begin{equation}
        a^\dagger_P
        a^\dagger_R
        a_S
        a_Q
    =
        -a^\dagger_P
        a^\dagger_R
        a_Q
        a_S
    =
        -a^\dagger_R
        a^\dagger_P
        a_S
        a_Q
\end{equation}
to remove any ambiguity about the sign we assume \( P < R \) and \(S > Q \)
and call this a canonical excitation.

In `NECI` any excitation of rank \(n\) is represented by a \(2 \times n\) matrix,
where the first row contains all the source indices (particles which are to be annihilated)
and the second row all target indices (particles which are to be created).

A canonical excitation \( a^\dagger_P a^\dagger_R a_S a_Q, (P < R \land Q < S )  \)
in `NECI` is given by a matrix, where each row is ascendingly sorted.
\begin{equation}
        \begin{bmatrix}
            \texttt{src}_1 & \texttt{src}_2 \\
            \texttt{tgt}_1 & \texttt{tgt}_2
        \end{bmatrix}
    =
        \begin{bmatrix}
            Q & S \\
            P & R
        \end{bmatrix}
    \quad, \quad
        \texttt{src}_1 < \texttt{src}_2, \quad \texttt{tgt}_1 < \texttt{tgt}_2
\end{equation}

This means that the operator, i.e. the matrix elements in `NECI`, is given by:
\begin{equation}
        \hat{g}
    =
        \sum_{PQRS}
        a^\dagger_P
        a^\dagger_R
        a_S
        a_Q
        g_{PQRS}
    =
        \sum_{PQRS}
        a^\dagger_P
        a^\dagger_R
        a_S
        a_Q
        U_{PRQS}
    =
        \sum
        \texttt{tgt}_1^\dagger \texttt{tgt}_2^\dagger \texttt{src}_2 \texttt{src}_1 U_{ \texttt{tgt}_1 \texttt{tgt}_2 \texttt{src}_2 \texttt{src}_1 }
\end{equation}
