module switch PrgEnv-pgi PrgEnv-gnu
module load scorep
module load vampir
export SCOREP_METRIC_PAPI=PAPI_L1_TCA,PAPI_L1_TCM,PAPI_TOT_INS
export SCOREP_ENABLE_TRACING=true
export SCOREP_TOTAL_MEMORY=1G
