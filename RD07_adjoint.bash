#!/bin/bash

#!/bin/bash

#====================
# QSUB
#====================

#$ -l gpu=1
#$ -l h_rt=10:50:0
#$ -l tmem=3G
#$ -N RTsolver
#$ -S /bin/bash

#$ -o RTsolver.txt
#$ -j y

#================================================================================
# Real Data 07 : finger
# 3D domain. 
# Compute the forward signal for sensors placed in the boundary of the cube
#================================================================================
if [ "$HOSTNAME" = "miller.local" ] || [ "$HOSTNAME" = "armstrong.local" ]; then
    export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build_miller/bin:$PATH"
else
    export PATH="/home/frullan/HighFreqCode/HighFreq_3DRT/Build/bin:$PATH"
fi

export EXAMPLE="RD07_finger3_doubleRes_randomSubsample/"
# Output folder
if [ "$HOSTNAME" = "maryam.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/ExperimentalData/"
elif [ "$HOSTNAME" = "ember.cs.ucl.ac.uk" ]; then
    export HOST_FOLDER="/scratch0/NOT_BACKED_UP/frullan/ExperimentalData/"
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
export INITIAL_PRESSURE="pixelPressure_0.dat"
export SENSORS="sensors_randomSubsampled_1800.dat" 
export SUBSAMPLE_INDICES="subsampleIndices_1800sensors.dat"
export FORWARD_SIGNAL="forwardSignal_random_1800sensors_390timesteps.dat"
export STDOUT="stdout-adjoint.txt"

# Mode
MODE='-a'
export GPU_INDEX=1
# Generate dimensions file
Nx=80  dx=0.000053
Ny=240 dy=0.000053
Nz=240 dz=0.000053
cat > $INPUT_FOLDER$DIMENSIONS <<EOF
$Nx $Ny $Nz
$dx $dy $dz
EOF

#==============================
# SENSORS
#==============================
nRaysPhi=1024 
nRaysTheta=1024
dt=1.6667e-8
tMax=6.49e-6
# Generate sensor file
cat > $INPUT_FOLDER$SENSORS <<EOF
EOF
# Step and tMax
echo "$dt $tMax 0 0 0 0 0 0 0" >> $INPUT_FOLDER$SENSORS
while read -r INDEX; 
    do 
    # Extract indices
    INDEX_Y=$(echo $INDEX | cut -d' ' -f 2)
    INDEX_Z=$(echo $INDEX | cut -d' ' -f 3)
    # Extract positions
    yPos=$(echo "scale=6;(2*($INDEX_Y-1)*$dy)" | bc)
    zPos=$(echo "scale=6;(2*($INDEX_Z-1)*$dz)" | bc)
    echo "0 $yPos $zPos $nRaysPhi $nRaysTheta -1.57 1.57 0.04 3.1" >> $INPUT_FOLDER$SENSORS;
done < $INPUT_FOLDER$SUBSAMPLE_INDICES


#==============================
# RUN ADJOINT
#==============================
RTsolver_GPU $MODE $INPUT_FOLDER$DIMENSIONS $INPUT_FOLDER$SOUND_SPEED $INPUT_FOLDER$INITIAL_PRESSURE \
             $INPUT_FOLDER$SENSORS $OUTPUT_FOLDER $INPUT_FOLDER$FORWARD_SIGNAL
