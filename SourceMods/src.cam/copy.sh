cp  ~/src/cam-opt/src/dynamics/se/dycore/control_mod.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' control_mod.F90
cp  ~/src/cam-opt/src/dynamics/se/dyn_grid.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' dyn_grid.F90
cp  ~/src/cam-opt/src/dynamics/se/dp_coupling.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' dp_coupling.F90
cp  ~/src/cam-opt/src/dynamics/se/dyn_comp.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' dyn_comp.F90
cp  ~/src/cam-opt/src/dynamics/se/dycore/fvm_mod.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' fvm_mod.F90
cp  ~/src/cam-opt/src/dynamics/se/dycore/prim_advection_mod.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' prim_advection_mod.F90
cp  ~/src/cam-opt/src/dynamics/se/dycore/global_norms_mod.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' prim_advection_mod.F90
cp  ~/src/cam-opt/src/dynamics/se/dycore/prim_driver_mod.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' prim_driver_mod.F90
cp  ~/src/cam-opt/src/dynamics/se/dycore/fvm_mapping.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' fvm_mapping.F90
cp  ~/src/cam-opt/src/dynamics/se/dycore/prim_advance_mod.F90 .
sed -i 's/!#define old_cam/#define old_cam/g' prim_advance_mod.F90
cp  ~/src/cam-opt/src/physics/cam/vertical_diffusion.F90 .

