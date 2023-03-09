#include <config.h>

void chanfast4(config_t *sim_config){
  validateConfig(sim_config);
  
  initMpiPlan();
  initDecompPlan();
  initFftwPlan();
  initTransposePlan();
  initGhostzPlan();
  initIo();
  initHelmholtzSolver();

  initGrid();

  allocateFlowMem();
  allocateStatsMem();

  initCfr();

  initFlow();

  for (; it < nt; it++){
    getCfl();
    for (step = 0; step < 3; step++){
      if (step == 0){
        resetRkAccumulation();
      }
      memSwap();

      ghostzPeriodic();
      applyBcs();

      solveFlowC();
      writeStatsC();

      solveFlowU();
      writeStatsU();
      
      solveFlowV();
      writeStatsV();

      solveFlowW();

      fractionalStep();
    }

    if ((it + 1) % istat == 0) {
      flushIo();
      writeFlowStatistics();
    }
    if ((it + 1) % iraw == 0) {
      writeFlowSnapshot();
    }
  }

  freeCfr();
  freeStatsMem();
  freeFlowMem();

  freeGrid();

  freeHelmholtzSolver();
  initIo();
  initGhostzPlan();
  initTransposePlan();
  initFftwPlan();
  initDecompPlan();
  initMpiPlan();

}