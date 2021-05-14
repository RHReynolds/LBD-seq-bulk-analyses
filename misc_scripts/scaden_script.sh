#!/bin/sh

#----Functions---------------------
function usage()
{
    echo "Script to run through scaden steps from creating bulk data to prediction"
    echo ""
    echo "-h --help"
    echo "--sample_dir=SAMPLE_DIR  Path to single-nuclear samples for bulk simulation"
    echo "--n_cells=CELLS  Number of cells to use in generation of bulk samples"
    echo "--n_simulated_samples=SIM_SAMPLES  Number of bulk simulated samples to generate per subject"
    echo "--n_training_steps=STEPS  Number of steps to train model with"
    echo "--pred_path=PRED_PATH  Path to prediction file"
    echo ""
}

#----Main--------------------------
while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        --sample_dir)
            SAMPLE_DIR=$VALUE
            ;;
        --n_cells)
            CELLS=$VALUE
            ;;
        --n_simulated_samples)
            SIM_SAMPLES=$VALUE
            ;;
        --n_training_steps)
            STEPS=$VALUE
            ;;
        --pred_path)
            PRED_PATH=$VALUE
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done


echo "sample_dir is $SAMPLE_DIR";
echo "n_cells is $CELLS";
echo "n_simulated_samples is $SIM_SAMPLES";
echo "n_training_steps is $STEPS";
echo "pred_path is $PRED_PATH";

# Navigate to where script should operate i.e. one directory path above where the sample data is stored
cd $SAMPLE_DIR
cd ..

# Make necessary directories
mkdir ./cells.${CELLS}_samples.${SIM_SAMPLES}/
mkdir ./cells.${CELLS}_samples.${SIM_SAMPLES}/samples_simulated/
mkdir ./cells.${CELLS}_samples.${SIM_SAMPLES}/scaden_model_${STEPS}steps/

# bulk simulation
python /opt/conda/pkgs/scaden-0.9.0-py_0/site-packages/scaden/preprocessing/bulk_simulation.py --cells $CELLS --samples $SIM_SAMPLES --data ${SAMPLE_DIR} --out ./cells.${CELLS}_samples.${SIM_SAMPLES}/samples_simulated/

# Create h5ad file
python /opt/conda/pkgs/scaden-0.9.0-py_0/site-packages/scaden/preprocessing/create_h5ad_file.py --data ./cells.${CELLS}_samples.${SIM_SAMPLES}/samples_simulated/ --out ./cells.${CELLS}_samples.${SIM_SAMPLES}/pdseq_bulksim.h5ad

# pre-processing of data
scaden process ./cells.${CELLS}_samples.${SIM_SAMPLES}/pdseq_bulksim.h5ad $PRED_PATH --processed_path ./cells.${CELLS}_samples.${SIM_SAMPLES}/processed_pdseq_bulksim.h5ad

# model training
scaden train  ./cells.${CELLS}_samples.${SIM_SAMPLES}/processed_pdseq_bulksim.h5ad --model_dir ./cells.${CELLS}_samples.${SIM_SAMPLES}/scaden_model_${STEPS}steps --steps  $STEPS

# predicting
scaden predict $PRED_PATH --model_dir ./cells.${CELLS}_samples.${SIM_SAMPLES}/scaden_model_${STEPS}steps --outname ./cells.${CELLS}_samples.${SIM_SAMPLES}/scaden_predictions_cells.${CELLS}:samples.${SIM_SAMPLES}.txt

