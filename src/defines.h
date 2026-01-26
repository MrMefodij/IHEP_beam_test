//
// Created by amefodev on 30.10.2024.
//

#ifndef CUBESLY_DEFINES_H
#define CUBESLY_DEFINES_H

#define AS_KV(x) #x << ": " << x
#define CHANNELS_MAP "../connection/connection_map.txt"
#define CUBES_MAP "../connection/cubes_vs_channels.txt"
#define CHANNELS_VS_CUBES "../connection/channels_vs_cubes.txt"
#define CHANNELS_NUMBER 50
#define CUBES_COUNT 27
#define LAYERS_COUNT 3
#define GLOBAL_CUBE 999 //3x3x3
#define CUBES_BOUNDARIES_MAP "../connection/cubes_boundaries.txt"
#define CUBES_SLICE_X 116
#define CUBES_SLICE_Y 116
#define CUBES_BINNING_X CUBES_SLICE_X, -33, 25
#define CUBES_BINNING_Y CUBES_SLICE_Y, -25, 33
#define CUBES_BINNING CUBES_BINNING_X, CUBES_BINNING_Y
#define FIBER_BINNING 600, 0, 300
#define GAUS_FIBER_BINNING 200, 0, 100
#define GAUS_BINNING 85, 0, 170
#define GAUS_BINNING_ANGLE 100, 20, 340
#define XTALK_BINNING 2*100*2.5, 0, 1
#define CUBES_POSITION 200
#define FIBRE_NOIZE 5000, 0, 200
#define CHANNELS_FLAG "../calibration/channel_flag.txt"
#define CUBES_NORMALIZATION "../calibration/cubes_normalization.txt"
#define CUBES_NORMALIZATION_ANGLE "../calibration/cubes_normalization_angle.txt"
#define FIBERS_NORMALIZATION "../calibration/fibers_normalization.txt"
#define FIBERS_NORMALIZATION_ANGLE "../calibration/fibers_normalization_angle.txt"
#define NOISE_SHIFT 11

#define CHARGE_FIT_CONSTANT 1.01275
#define DEFAULT_BINNING 21
#define SINGLE_CUBE_BINNING_1D 24, 0, 12
#define SINGLE_CUBE_BINNING SINGLE_CUBE_BINNING_1D, SINGLE_CUBE_BINNING_1D
#define AMPLITUDE_OVERFLOW 3000
#define XTALK_BINNING_MEAN 1000, 0 , 1
#define AXIS_RANGE_XTALK_2D_Z 0.02, 0.08
#define XTALK_RANGE_1D 0, 0.1
#endif //CUBESLY_DEFINES_H
