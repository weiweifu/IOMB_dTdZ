#!/bin/sh

# ~/ilamb/bin/ilamb-run --config my.cfg --model_root $ILAMB_ROOT/MODELS/ --regions global  --models CanESM5 --clean 

# ~/ILAMB/bin/ilamb-run --config my.cfg --model_root $ILAMB_ROOT/MODELS/ --regions global  

~/ilamb/bin/ilamb-run --config m0.cfg --model_root $ILAMB_ROOT/MODELS/ --regions global --build_dir test03 --models CNRM-ESM2 UKESM --clean 



# ~/ILAMB/bin/ilamb-run --config new.cfg --model_root $ILAMB_ROOT/MODELS/ --regions global --build_dir testm0 --models  MPI-ESM1-2-HR CNRM-ESM2 --clean 

# ~/ilamb/bin/ilamb-run --config new3.cfg --model_root $ILAMB_ROOT/MODELS/ --regions global --build_dir IOMB_CMIP6 --model_setup model_order.txt 
# --models CNRM-ESM2 UKESM --clean 
