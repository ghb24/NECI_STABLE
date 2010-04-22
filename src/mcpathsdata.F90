module mcpathsdata
       TYPE EGP
         integer, pointer :: p(:) !memory store
         integer          :: l    !length of memory store
         integer          :: v    !vertex mem store was from
       end type
end module mcpathsdata
