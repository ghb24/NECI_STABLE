! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module mcpathsdata
       TYPE EGP
         integer, pointer :: p(:) !memory store
         integer          :: l    !length of memory store
         integer          :: v    !vertex mem store was from
       end type
end module mcpathsdata
