How to use CubesLY
1) Each run should have calibration results with name "testbench_#_Calibration_A" & "testbench_#_Calibration_I" in folder "./calibration/A/" and "./calibration/I/";
2) If You don't have that you can do that with for runs testbench_601 to testbench_609
   `for n in {1..9}; do;./mppcCalibration ~/work/data/PNPI_2024/new_unpack/testbench_60${n}.root; done;`
3) To do conversion from PNPI structure to hitSuperstructure please do for runs testbench_601 to testbench_609:
    `for n in {1..9}; do;./eventStructure ~/work/data/PNPI_2024/new_unpack/testbench_60${n}.root; done;`
     OutputFile in data folder : `testbench_60${n}_HitStructure.root` & `testbench_#_Calibration_A.txt` & `testbench_#_Calibration_I.txt`
4) To calculate cubes boundaries please use:
    `./cubesLocation /Users/amefodev/work/data/PNPI_2024/new_unpack/testbench_605.root ` & `../connection/channels_vs_cubes.txt`
    OutputFile in data folder : `All_runs_cubes_position.root`
5) To make plots for LY you need all calibration file and hitStructure prepared, that macros use all `*_HitStructure.root` files in folder,
    `./cubesResponse /Users/amefodev/work/data/PNPI_2024/new_unpack/testbench_601_HitStructure.root `
6) To run xTalk macros you need all calibration file and hitStructure prepared, that macros use all `*_HitStructure.root` files in folder,
    `./xTalk /Users/amefodev/work/data/PNPI_2024/new_unpack/testbench_601_HitStructure.root `