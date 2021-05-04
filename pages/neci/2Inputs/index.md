title: Calculation Inputs 
--- 

## Calculation inputs

The NECI executable takes one input argument, which is the name of an
input file containing the instructions for carrying out the calculation.
The input file is organized in blocks, with each block being started and
terminated by a dedicated keyword. Each block can contain a number of
keywords to specify options. Here, a list of the blocks and their
respective keywords is given.

The first line of the input is always `title`, the last line is always
`end`.

Some keywords are mandatory, those are marked in <span
style="color: red">**red**</span> and are given at the beginning of the
description of each paragraph. Then come recommended options, marked in
<span style="color: blue">**blue**</span>, followed by further options
given in black.

Keywords which are purely for debugging purposes and only interesting
for developers are markes as <span style="color: green">**green**</span>.

- [ ] TODO I am not sure if this should just one page or several 