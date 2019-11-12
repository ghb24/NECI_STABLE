#:def _get_decl(type, kind=None)
$:  '{}{}'.format(type, '('+kind+')' if kind else '')
#:enddef _get_decl

#:def _get_name(type, kind=None)
$:  '{}{}'.format(type, '_'+kind if kind else '')
#:enddef _get_name

#:def get_sort_name(type, kind=None)
$:  'sort_'+_get_name(type, kind)
#:enddef get_sort_name

#:def create_sort(type, kind=None)
    function ${get_sort_name(type, kind)}$(x) result(res)
        ${_get_decl(type, kind)}$, intent(in) :: x(:)
        ${_get_decl(type, kind)}$ :: res

        res = x(1)
    end function ${get_sort_name(type, kind)}$
#:enddef create_sort

