#!/bin/bash

#====================
# QSUB
#====================

# #$ -P gpu
#$ -l gpu=1
#$ -l h_rt=200:00:00
#$ -l tmem=3G
#$ -N fistaRD10_tau1e1
#$ -wd /home/frullan/HighFreqCode/ExperimentalData/RD10_finger2_doubleRes_subsampled
#$ -S /bin/bash

# -o RTiter.txt
#$ -j y

#================================================================================
# EXAMPLE 86 ITERATIVE
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
#The code you want to run now goes here.
if [ "$HOSTNAME" = "miller.local" ] || [ "$HOSTNAME" = "armstrong.local" ]; then
    export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build_miller/bin:$PATH"
else
    export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
fi

export EXAMPLE="RD10_finger2_doubleRes_subsampled/"

# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/"
elif [ "$HOSTNAME" = "hannover" ]; then
    export HOST_FOLDER="/home/wontek/sharedWK/ExperimentalData/"
else
    export HOST_FOLDER="/home/frullan/HighFreqCode/ExperimentalData/"
fi
export EXAMPLE_FOLDER=$HOST_FOLDER$EXAMPLE
export INPUT_FOLDER=$EXAMPLE_FOLDER"input_data/"
export OUTPUT_FOLDER=$EXAMPLE_FOLDER"output_data/"
cd $EXAMPLE_FOLDER

# Assign files
export DIMENSIONS="dimensions.dat"
export SOUND_SPEED="sound_speed.dat"
export INITIAL_PRESSURE="initial_pressure_veins_80x240x240.dat"
export SENSORS="sensors_subsampled_3600.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_3600sensors_400timesteps.dat"
export PIXEL_PRESSURE="pixelPressure_0.dat"

# Machine
echo $HOSTNAME
# Choose GPU
export GPU_INDEX=0
# Choose mode
export MODE='-F'

# Parameters
SIGMA=5e-2
TAU=1e1
THETA=1    
LAMBDA=5e-5
BATCH_SIZE=100
NITER=200

#=======   GRADIENT DESCENT  ========================================
if [ "$MODE" = "-G" ]; then
    echo "=================== GRADIENT DESCENT ===================="
    export STDOUT="stdout_GD_tau"$TAU"_lambda"$LAMBDA$"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER
# > $OUTPUT_FOLDER$STDOUT
#=======   STOCHASTIC GRADIENT DESCENT
elif [ "$MODE" = "-g" ]; then
    echo "=================== STOCHASTIC GRADIENT DESCENT ===================="
    export STDOUT="stdout_S-GD_tau"$TAU"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $BATCH_SIZE $NITER > $OUTPUT_FOLDER$STDOUT
#=======   FISTA
elif [ "$MODE" = "-F" ]; then
    echo "=================== FISTA ===================="
    export STDOUT="stdout_FISTA_tau"$TAU"_lambda"$LAMBDA"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#=======   PRIMAL DUAL HYBRID GRADIENT
elif [ "$MODE" = "-P" ]; then
    echo "=================== PDHG ===================="
    export STDOUT="stdout_PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#=======   STOCHASTIC PRIMAL DUAL HYBRID GRADIENT
elif [ "$MODE" = "-p" ]; then
    echo "=================== S-PDHG ===================="
    export STDOUT="stdout_S-PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $BATCH_SIZE $NITER > $OUTPUT_FOLDER$STDOUT
#=======   SINGLE FORWARD ADJOINT
elif [ "$MODE" = "-r" ]; then
    echo "============  SINGLE FORWARD ADJOINT  ============"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE
else
    echo "Non supported mode"
fi

    
