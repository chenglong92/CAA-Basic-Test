objects=DRP7_sub.f90 f_sub.f90 LDDRK_sub.f90 Lax_Wendroff_sub.f90 \
	One_dimensional_migration_main.f90
DRP.exe: $(objects)
	gfortran $(objects) -o DRP.exe
