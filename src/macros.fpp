#:def unused_var(*args)
#ifdef WARNING_WORKAROUND_
    #:for arg in args
        associate(${arg}$ => ${arg}$); end associate
    #:endfor
#endif
#:enddef unused_var
