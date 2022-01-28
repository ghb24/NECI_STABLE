!Original author: Anthony Stone, since modified by ajwt and ghb.

MODULE input_neci

    use constants, only: sp, dp, int64, stdout
    use util_mod, only: stop_all, get_free_unit

    IMPLICIT NONE

    integer, parameter :: max_len = 300
    CHARACTER(LEN=1900) :: char = ""
    logical, parameter :: skipbl = .true., clear = .true.
    LOGICAL :: echo = .false., more
    INTEGER :: item = 0, nitems = 0, loc(0 : max_len) = 0, end(max_len) = 0, &
                     line = 0, nerror = 0, ir = 5, last = 0, unit(0:10)

    CHARACTER(LEN=26), PARAMETER :: &
        upper_case = "ABCDEFGHIJKLMNOPQRSTUVWXYZ", &
        lower_case = "abcdefghijklmnopqrstuvwxyz"
    CHARACTER(len=1), PARAMETER :: space = " ", bra = "(", ket = ")", &
                            comma = ",", squote = "'", dquote = '"', tab = achar(9), &
                            plus = "+", minus = "-", dot = "."

    CHARACTER(len=*), parameter :: concat = "+++"
    INTEGER, parameter :: lc = len(concat)
    CHARACTER(LEN=455) :: file = ""

    external :: neci_getarg

    interface readf
        module procedure read_single
        module procedure read_double
    end interface

    save
    PRIVATE
    public :: item, nitems, read_line, input_options, readu, geti, ir, char, &
        readi, readf, readl, reada, reread, report, getilong, getf, getrange
!  AJWT - added , ir to above
!  Free-format input routines

!     CALL READ_LINE(eof[,inunit])
!  Read next input record from unit IR into buffer, and analyse for data
!  items, null items (enclosed by commas), and comment (enclosed in
!  parentheses, which may be nested, and occurring where a new item may
!  occur). Data items are terminated by space or comma, unless enclosed
!  in single or double quotes.
!  If the optional argument inunit is provided, the record is read from
!  there instead.
!  If a line ends with the concatenation sequence (default "+++")
!  possibly followed by spaces, the sequence and any following spaces
!  are removed and the next line concatenated in its place.
!  This is repeated if necessary until a line is read that does not end
!  with the concatenation sequence.
!  The logical variable eof is set TRUE if end of file is encountered,
!  otherwise FALSE.
!  The public module variable ITEM is set to zero and NITEMS to the number
!  of items in the record.

!     CALL PARSE(string)
!  Parse the string provided in the same way as the read_line routine
!  (which uses this routine itself) and leave the details in the buffer
!  as for a line read from the input. Input directives are not
!  interpreted.

!     CALL READx(V)
!  Read an item of type x from the buffer into variable V:
!     CALL READF   single or double precision, depending on the type of V
!     CALL READI   integer
!     CALL READILONG integer(int64)
!     CALL READA   character string
!     CALL READU   character string, uppercased
!     CALL READL   character string, lowercased
!  All of these return null or zero if there are no more items on the
!  current line. Otherwise ITEM is incremented so that it gives the
!  number of the item that has just been read.
!     CALL READF(V,factor)
!  If the optional second argument is supplied to READF, the variable V
!  is divided by it after it has been read. This is convenient for converting
!  from external to internal units.
!     CALL REREAD(K)
!  K>0  Prepare to read item K
!  K<0  Go back |K| items
!  K=0  Same as K=-1, i.e. reread last item.

!     CALL GETx
!  Same as the corresponding READx, but a new record is read if there are
!  no more items in the current one.

!  ITEM    is the number of the last item read from the buffer

!  NITEMS  is the number of items in the buffer

!  char(I:I) is the Ith character in the buffer

!  LOC(I)  is the starting position of the Ith item
!  END(I)  is the position of the last character of the Ith item

!  LINE    is the number of lines (records) read

!  If SKIPBL is set to TRUE, lines containing no items other than
!  comment are skipped automatically, so that INPUT will always return
!  a line with at least one item (possibly null).  If SKIPBL is FALSE
!  (default) no skipping occurs and the next data line is returned
!  regardless of what it contains.

!  If CLEAR is set to TRUE (default) then null items will be returned
!  as zero or blank. If an attempt is made to read more than
!  NITEMS items from a line, the items are treated as null. If
!  CLEAR is FALSE, a variable into which a null item is read is
!  left unaltered.

!  NERROR specifies the treatment of errors found while reading numbers:
!    0  Hard error - print message and stop (default).
!    1  Soft error - print message and return zero.
!    2  Soft error - no message, return zero, set NERROR to -1.
!  If NERROR is set to 2 it is important to test and reset it after reading
!  a number.

!  IR is the input stream from which the data are to be read (default 5).

!  If ECHO is TRUE, the input line will be reflected to standard output.

!  LAST gives the position of the last non-blank character in the line.

!     CALL REPORT(message[,REFLECT])
!  Error-message routine. Prints the error message, followed by the
!  current input buffer if REFLECT is present and TRUE, and stops.

CONTAINS

    SUBROUTINE read_line(eof, inunit)

!  Read next input record from unit IR into buffer, and analyse for data
!  items, null items (enclosed by commas), and comment (enclosed in
!  parentheses, which may be nested, and occurring where spaces may
!  occur).
!  If the optional argument inunit is specified, read a line from that
!  unit instead of IR.

!  Stream-switching commands may occur in the data:
!      #include file-name
!         switches input to be from the file specified;
!      #revert
!         (or end-of-file) reverts to the previous file.
!  However it does not make sense for stream-switching commands to
!  be used when the input unit is specified explicitly. If they are given
!  in this case, they apply to the default input stream.
!    Also
!      #concat "string"
!         sets the concatenation string; e.g.
!         #concat "\"
!      #width <n>
!         Input line width limited to n characters. Default is 128; may
!         need to be reduced to 80 for some IBM computers. Maximum is 255.

!  CONCAT is the line concatenation string: if it is found as the last
!  non-blank character string on a line, the following line is read and
!  appended to the current line, overwriting the concatenation string.
!  This procedure is repeated if necessary until a line is found that does
!  not end with the concatenation string.

!  The concatenation string may be modified by changing its declaration
!  above. If it is null, no concatentation occurs. The concatenation string
!  may also be changed by the #concat directive in the data file.

        LOGICAL, INTENT(OUT) :: eof
        INTEGER, INTENT(IN), OPTIONAL :: inunit

        character(*), parameter :: this_routine = 'read_line'

        CHARACTER(LEN=455) :: w, f
        CHARACTER :: term

        INTEGER, SAVE :: lrecl = 466
        INTEGER :: in, fail, i, k, l, m, ierr

        eof = .false.
        if (present(inunit)) then
            in = inunit
        else
            in = ir
        end if

        char = ""
        lines: do
            more = .true.
            m = 1
            do while (more)
                last = m + lrecl - 1
                line = line + 1
                read(in, "(a)", iostat=ierr) char(m:last)
                if (ierr > 0) then
                    call stop_all(this_routine, 'Error during read')
                else if (ierr < 0) then
                    ! EOF reached
                    eof = .true.
                    exit
                else
    !  Find last non-blank character
                    last = verify(char, space//tab, back=.true.)
                    if (echo) write(stdout, "(a)") char(m:last)
    !  Look for concatenation string
                    if (lc > 0 .and. last >= lc) then
                        more = (char(last - lc + 1:last) == concat)
                        if (more) then
                            m = last - lc + 1
                        end if
                    else
                        more = .false.
                    end if
                end if
            end do  ! while (more)

            is_EOF_reached: if (.not. eof) then

    !  Replace tab by single space
                do while (index(char, tab) > 0)
                    L = index(char, tab)
                    char(L:L) = space
                end do

    !  Logical line assembled. First look for input directives
                L = 1
                do while (char(L:L) == space .and. L < last)
                    L = L + 1
                end do
                if (char(L:L) == "#") then
                    M = L
                    do while (char(M:M) /= space .and. M <= last)
                        M = M + 1
                    end do
                    w = char(L:M - 1)
                    call upcase(w)
                    if (M > last) then
                        f = " "
                    else
                        do while (char(M:M) == space)
                            M = M + 1
                        end do
                        if (char(M:M) == squote .or. char(M:M) == dquote) then
                            term = char(M:M)
                            M = M + 1
                        else
                            term = space
                        end if
                        L = M
                        do while (char(M:M) /= term .and. M <= last)
                            M = M + 1
                        end do
                        f = char(L:M - 1)
                    end if
                    cycle lines
                end if

                call parse

    !  Blank except for comment?
                if (nitems == 0 .and. skipbl) then
                    cycle lines   !  Read another line
                else
                    exit lines    !  Finished
                end if
            else
    !  End of file
                if (more .and. m > 1) then
                    write(stdout, "(a)") "Apparently concatenating at end-of-file"
                    call report("Unexpected end of data file", .true.)
                else
                    !  End of input
                    char(1:last) = ""
                    item = 0
                    nitems = -1
                    eof = .true.
                    return
                end if
            end if is_EOF_reached
        end do lines
    END SUBROUTINE read_line

!-----------------------------------------------------------------------

    SUBROUTINE parse(string)

        CHARACTER(LEN=*), OPTIONAL :: string

        INTEGER :: L, state, nest
        LOGICAL :: tcomma
        CHARACTER :: term, c

        if (present(string)) then
            char = string
            last = len_trim(string)
        end if

!  Analyse input
!  State numbers:
!  0  Looking for item
!  1  Reading quoted string
!  2  Reading unquoted item
!  3  Reading comment
!  4  Expecting space or comma after quoted string

        state = 0
        item = 0
        nitems = 0
        L = 0            ! Position in input buffer
        tcomma = .true.  !  True if last item was terminated by comma
        !  and also at start of buffer

        chars: do

!  Read next character
            L = L + 1
            if (L > last) then
!  End of line
                if (nitems > 0) then
                    select case (state)
                    case (1)
                        call report("Closing quote missing", .true.)
                    case (2)
                        end(nitems) = L - 1
                    end select
                end if
                return
            end if

            c = char(L:L)
            select case (state)
            case (0)                ! Looking for next item
                select case (c)
                case (space, tab)      ! Keep on looking
                    !           cycle chars
                case (bra)            ! Start of comment
                    nest = 1
                    state = 3
                    !           cycle chars
                case (squote, dquote)  ! Start of quoted string
                    nitems = nitems + 1
                    loc(nitems) = L
                    term = c
                    state = 1
                case (comma)
                    if (tcomma) then   ! Null item between commas
                        nitems = nitems + 1
                        loc(nitems) = 0
                    end if
                    tcomma = .true.
                case default         ! Start of unquoted item
                    nitems = nitems + 1
                    loc(nitems) = L
                    state = 2
                end select

            case (1)                ! Reading through quoted string
                if (c == term) then ! Closing quote found
                    end(nitems) = L
                    state = 4
                    !         else
                    !           cycle chars
                end if

            case (2)                ! Reading through unquoted item
                select case (c)
                case (space, tab)      ! Terminator
                    end(nitems) = L - 1
                    state = 0
                    tcomma = .false.
!   case(bra)            ! Start of comment -- treat as space
!     end(nitems)=L-1    ! This code allows for parenthesised comments
!     tcomma=.false.     ! to be embedded in unquoted strings. This is
!     state=3            ! not a good idea. Such comments now have to occur
!     nest=1             ! only where a new item might begin.
                case (comma)          ! Comma-terminated
                    end(nitems) = L - 1
                    state = 0
                    tcomma = .true.
                    !         case default
                    !           cycle chars
                end select

            case (3)                ! Reading through comment
                select case (c)
                case (bra)            ! Nested parenthesis
                    nest = nest + 1
                case (ket)
                    nest = nest - 1
                    if (nest == 0) then  ! End of comment
                        state = 0          ! Space or comma not required after comment
                    end if
                    !         case default
                    !           cycle chars
                end select

            case (4)                ! Expecting space or comma
                select case (c)
                case (space, tab)
                    tcomma = .false.
                    state = 0
                case (comma)
                    tcomma = .true.
                    state = 0
                case (bra)            ! Start of comment -- treat as space
                    tcomma = .false.
                    state = 3
                    nest = 1
                case default
                    call report("Space or comma needed after quoted string", .true.)
                end select

            end select

        end do chars

    END SUBROUTINE parse

!-----------------------------------------------------------------------

    SUBROUTINE input_options(echo_lines)

        LOGICAL, intent(in) :: echo_lines

        echo = echo_lines
    END SUBROUTINE input_options

!-----------------------------------------------------------------------

    SUBROUTINE reada(m)

!  Copy characters from the next item into the character variable M.
!  If the first character is a single or double quote, the string is
!  terminated by a matching quote and the quotes are removed.

        CHARACTER(LEN=*), INTENT(INOUT) :: m
        INTEGER :: l

        if (clear) m = ""
!  If there are no more items on the line, M is unchanged
        if (item >= nitems) return

        item = item + 1
!  Null item?
        if (loc(item) == 0) return

        l = loc(item)
        if (char(l:l) == squote .or. char(l:l) == dquote) then
            m = char(l + 1:end(item) - 1)
        else
            m = char(l:end(item))
        end if

    END SUBROUTINE reada

!-----------------------------------------------------------------------

    SUBROUTINE read_double(A, factor)

!  Read the next item from the buffer as a real (double precision) number.
!  If the optional argument factor is present, the value read should be
!  divided by it. (External value = factor*internal value)

        REAL(KIND=dp), INTENT(INOUT) :: a
        REAL(KIND=dp), INTENT(IN), OPTIONAL :: factor
        CHARACTER(LEN=90) :: string

        if (clear) a = 0d0

!  If there are no more items on the line, I is unchanged
        if (item >= nitems) return

        string = ""
        call reada(string)
!  If the item is null, I is unchanged
        if (string == "") return
        read(unit=string, fmt=*, err=99) a
        if (present(factor)) then
            a = a / factor
        end if
        return

99      a = 0d0
        select case (nerror)
        case (-1, 0)
            call report("Error while reading real number", .true.)
        case (1)
            write(stdout, "(2a)") "Error while reading real number. Input is ", trim(string)
        case (2)
            nerror = -1
        end select

    END SUBROUTINE read_double

!-----------------------------------------------------------------------

    SUBROUTINE read_single(A, factor)

!  Read the next item from the buffer as a real (double precision) number.
!  If the optional argument factor is present, the value read should be
!  divided by it. (External value = factor*internal value)

        REAL(kind=sp), INTENT(INOUT) :: a
        REAL(kind=sp), INTENT(IN), OPTIONAL :: factor

        REAL(kind=dp) :: aa

        if (present(factor)) then
            call read_double(aa, real(factor, dp))
        else
            call read_double(aa)
        end if
        a = real(aa, sp)

    END SUBROUTINE read_single

!-----------------------------------------------------------------------

    SUBROUTINE readi(I)
!  Read an integer from the current record

        INTEGER, INTENT(INOUT) :: i

        CHARACTER(LEN=90) :: string

        if (clear) i = 0

!  If there are no more items on the line, I is unchanged
        if (item >= nitems) return

        string = ""
        call reada(string)
!  If the item is null, I is unchanged
        if (string == "") return
        read(unit=string, fmt=*, err=99) i
        return

99      i = 0
        select case (nerror)
        case (-1, 0)
            call report("Error while reading integer", .true.)
        case (1)
            write(stdout, "(2a)") "Error while reading integer. Input is ", trim(string)
        case (2)
            nerror = -1
        end select

    END SUBROUTINE readi

!---------------------------------------------------

    SUBROUTINE readiLong(I)
!  Read a long integer from the current record

        integer(int64), INTENT(INOUT) :: i

        CHARACTER(LEN=90) :: string

        if (clear) i = 0

!  If there are no more items on the line, I is unchanged
        if (item >= nitems) return

        string = ""
        call reada(string)
!  If the item is null, I is unchanged
        if (string == "") return
        read(unit=string, fmt=*, err=99) i
        return

99      i = 0
        select case (nerror)
        case (-1, 0)
            call report("Error while reading long integer", .true.)
        case (1)
            write(stdout, "(2a)") "Error while reading long integer. Input is ", trim(string)
        case (2)
            nerror = -1
        end select

    END SUBROUTINE readiLong

!-----------------------------------------------------------------------

    SUBROUTINE readu(m)
        CHARACTER(LEN=*) m

        call reada(m)
        call upcase(m)

    END SUBROUTINE readu

!-----------------------------------------------------------------------

    SUBROUTINE readl(m)
        CHARACTER(LEN=*) m

        call reada(m)
        call locase(m)

    END SUBROUTINE readl

!-----------------------------------------------------------------------

    SUBROUTINE getf(A, factor)
!  Read the next item as a double-precision number, reading new data
!  records if necessary.
!  If the optional argument factor is present, the value read should be
!  divided by it. (External value = factor*internal value)

        REAL(kind=dp), INTENT(INOUT) :: A
        REAL(kind=dp), INTENT(IN), OPTIONAL :: factor

        LOGICAL :: eof

        do
            if (item < nitems) then
                call readf(a, factor)
                exit
            else
                call read_line(eof)
                if (eof) then
                    call stop_all('getf', "End of file while attempting to read a number")
                end if
            end if
        end do

    END SUBROUTINE getf

!-----------------------------------------------------------------------

    SUBROUTINE geti(I)
!  Get an integer, reading new data records if necessary.
        INTEGER, INTENT(INOUT) :: i
        LOGICAL :: eof

        do
            if (item < nitems) then
                call readi(i)
                exit
            else
                call read_line(eof)
                if (eof) then
                    call stop_all('geti', "End of file while attempting to read a number")
                end if
            end if
        end do

    END SUBROUTINE geti

!-----------------------------------------------------------------------

    SUBROUTINE getiLong(I)
!  Get an integer, reading new data records if necessary.
        integer(int64), INTENT(INOUT) :: i
        LOGICAL :: eof

        do
            if (item < nitems) then
                call readiLong(i)
                exit
            else
                call read_line(eof)
                if (eof) then
                    call stop_all('getiLong', "End of file while attempting to read a number")
                end if
            end if
        end do

    END SUBROUTINE getiLong

!-----------------------------------------------------------------------


    SUBROUTINE reread(k)

        INTEGER, INTENT(IN) :: k
!  k>0  Reread from item k
!  k<0  Go back |k| items
!  k=0  Same as k=-1, i.e. reread last item.

        if (k < 0) then
            item = item + k
        else if (k == 0) then
            item = item - 1
        else
            item = k - 1
        end if
        if (item < 0) item = 0

    END SUBROUTINE reread

!-----------------------------------------------------------------------

    subroutine getRange(w, start, end)
        character(*), intent(inout) :: w
        integer, intent(out) :: start, end
        integer :: index

        w = adjustl(trim(w))
        index = scan(w, "-")
        read(w(1:index - 1), *) start
        read(w(index + 1:), *) end
    end subroutine getRange

!-----------------------------------------------------------------------

    SUBROUTINE upcase(word)
        CHARACTER(LEN=*), INTENT(INOUT) :: word
        INTEGER :: i, k

        do i = 1, len(word)
            k = index(lower_case, word(i:i))
            if (k /= 0) word(i:i) = upper_case(k:k)
        end do

    END SUBROUTINE upcase

!-----------------------------------------------------------------------

    SUBROUTINE locase(word)
        CHARACTER(LEN=*), INTENT(INOUT) :: word
        INTEGER :: i, k

        do i = 1, len(word)
            k = index(upper_case, word(i:i))
            if (k /= 0) word(i:i) = lower_case(k:k)
        end do

    END SUBROUTINE locase

!-----------------------------------------------------------------------

    SUBROUTINE report(c, reflect)

        CHARACTER(LEN=*), INTENT(IN) :: c
        LOGICAL, INTENT(IN), OPTIONAL :: reflect
        INTEGER :: i, i1, i2, l

        CHARACTER(LEN=3) s1, s2

        write(stdout, "(a)") c
        if (present(reflect)) then
            if (reflect) then
                l = loc(item)
                i2 = min(last, l + 20)
                i1 = max(1, i2 - 70)
                s1 = " "
                if (i1 > 1) s1 = "..."
                s2 = " "
                if (i2 < last) s2 = "..."
                write(stdout, "(a, I5)") "Input line ", line
                write(stdout, "(a3,1x,a,1x,a3)") s1, char(i1:i2), s2
                write(stdout, "(3x,80a1)") (" ", i=i1, l), "*"
            end if
        end if
        call stop_all("report", 'Input error')

    END SUBROUTINE report


END MODULE input_neci
