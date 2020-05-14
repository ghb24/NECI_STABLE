! A sensible way to ensure that the current case name is always set for
! any tests
!
! Replace: call test_fn(args) --> TEST1(test_fn, args)
! Replace: call test_fn --> TEST(test_fn)

#define TESTs1(name,subname,args) \
	call set_case_name(__FILE__ // ": " // "name" // " (" // trim(adjustl(subname)) // ")"); \
	call name(args); \
	call set_case_name('_not_set_')

#define TEST1(name,args) \
	call set_case_name(__FILE__ // ": " // "name"); \
	call name(args); \
	call set_case_name('_not_set_')

#define TEST(name) \
	call set_case_name(__FILE__ // ": " // "name"); \
	call name; \
	call set_case_name('_not_set_')

