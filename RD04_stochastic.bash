#!/bin/bash

#====================
# QSUB
#====================

#$ -P gpu
#$ -l gpu=1
#$ -l h_rt=20:00:00
#$ -l tmem=3G
#$ -N spdhg04_tau2e18_sigma1_lambda5e-5
#$ -wd /home/frullan/HighFreqCode/ExperimentalData/RD04_palm/
#$ -S /bin/bash

# -o RTiter.txt
#$ -j y


#================================================================================
# RD04 Real data
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================

export EXAMPLE="RD04_palm/"

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
export SENSORS="sensors_subsampled_8736.dat" 
export FORWARD_SIGNAL="forwardSignal_reference_norm_8736sensors.dat"
export PIXEL_PRESSURE="pixelPressure_0.dat"

# Choose GPU
export GPU_INDEX=1
# Choose mode
export MODE='-G'

#================================================================================
#=======   GRADIENT DESCENT
#================================================================================
if [ "$MODE" = "-G" ]; then
    echo "=================== GRADIENT DESCENT ===================="
    # Regularization parameters
    TAU=5e17
    LAMBDA=2e-3
    NITER=10
    # Output
    export STDOUT="stdout_GD_tau"$TAU"_lambda"$LAMBDA$"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   STOCHASTIC GRADIENT DESCENT
#================================================================================
elif [ "$MODE" = "-g" ]; then
    echo "=================== STOCHASTIC GRADIENT DESCENT ===================="
    # Regularization parameters
    TAU=2e18
    LAMBDA=1e-3
    BATCH_SIZE=800
    N_EPOCHS=10
    # Output
    export STDOUT="stdout_S-GD_tau"$TAU"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$N_EPOCHS".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$SENSORS \
                    $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $BATCH_SIZE $N_EPOCHS > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   FISTA
#================================================================================
elif [ "$MODE" = "-F" ]; then
    echo "=================== FISTA ===================="
    # Regularization parameters
    TAU=2e17
    LAMBDA=1e-2
    NITER=10
    # Output
    export STDOUT="stdout_FISTA_tau"$TAU"_lambda"$LAMBDA"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   STOCHASTIC FISTA
#================================================================================
elif [ "$MODE" = "-f" ]; then
    echo "=================== STOCHASTIC FISTA ===================="
    # Regularization parameters
    TAU=1e18   
    LAMBDA=1e-3
    BATCH_SIZE=90
    N_EPOCHS=10
    # Output
    export STDOUT="stdout_S-FISTA_tau"$TAU"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$N_EPOCHS".txt"
    RTiterative_GPU $1 $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $TAU $LAMBDA $BATCH_SIZE $N_EPOCHS> $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   PRIMAL DUAL HYBRID GRADIENT
#================================================================================
elif [ "$MODE" = "-P" ]; then
    echo "=================== PDHG ===================="
    # Regularization parameters
    SIGMA=5e-1
    TAU=5e17
    THETA=1      
    LAMBDA=1e-3  
    NITER=10
    # Output
    export STDOUT="stdout_PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_iter"$NITER".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $NITER > $OUTPUT_FOLDER$STDOUT
#================================================================================
#=======   STOCHASTIC PRIMAL DUAL HYBRID GRADIENT
#================================================================================
elif [ "$MODE" = "-p" ]; then
    echo "=================== S-PDHG ===================="
    # Regularization parameters
    SIGMA=1
    TAU=2e18
    THETA=1    
    LAMBDA=5e-5
    BATCH_SIZE=90
    N_EPOCHS=10
    # Output
    export STDOUT="stdout_S-PDHG_sigma"$SIGMA"_tau"$TAU"_theta"$THETA"_lambda"$LAMBDA"_batch"$BATCH_SIZE"_epochs"$N_EPOCHS".txt"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE $SIGMA $TAU $THETA $LAMBDA $BATCH_SIZE $N_EPOCHS > $OUTPUT_FOLDER$STDOUT
elif [ "$MODE" = "-r" ]; then
    echo "============  SINGLE FORWARD ADJOINT  ============"
    RTiterative_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED \
                    $INPUT_FOLDER$SENSORS $INPUT_FOLDER$FORWARD_SIGNAL $INPUT_FOLDER$PIXEL_PRESSURE
else
    echo "Non supported mode"
fi

    
