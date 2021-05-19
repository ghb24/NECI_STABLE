[^keygen]: `ssh-keygen` can also generate DSA keys. Some ssh clients and
servers will reject DSA keys longer than 1024 bits, and 1024 bits is
currently on the margin of being crackable. As such 2048 bit RSA keys
are preferred. Top secret this code is. Probably. Apart from the master
branch which hosted for all on github. And in molpro. And anyone that
wants it obviously.

[^powerpitzer]: V. Neufeld, A. Thom, J. Chem. Theory Comput.2019151127-140

[^initiator]: D. Cleland, G.H. Booth, A. Alavi, J. Chem. Phys. 132, 041103 (2010)

[^kpfciqmc]: N. S. Blunt, Ali Alavi, George H. Booth, Phys. Rev. Lett. 115,
050603

[^ctags]: Executing `ctags --version` should print either `Universal Ctags` or
`Exuberant Ctags` but not `ctags (GNU Emacs)`

[^gethelement]: As an exception, some old code makes use of `gethelement` or
`gethelement2`. These should not be used in new code. Some of the
initialisation code also uses `gethelement` as the prerequisites for
`get_helement` have not yet been met.

[^diagelem]: This general routine is not implemented at the time of writing, but
will added shortly.

[^double_excit]: Note that for double excitations, in principle there are two
alignments that work. The two new orbitals could be either way around
and the parity of these versions are inverted relative to each other. A
convention for the ordering of the new orbitals (in NECI they are
required to be numerically increasing) is required but arbitrary. The
overall simulation will give the same results either way.

[^test]: test if this is included.