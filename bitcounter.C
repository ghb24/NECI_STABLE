//This should be a quicker bit counter than the other one, but try to convert to fortran.
int NumberOfSetBits(int i)

{

    i = i - ((i >> 1) & 0x55555555);

    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);

    return ((i + (i >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;

}
