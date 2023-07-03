#include <check.h>
#include <logging.h>
#include <io.h>

#include <globals_sim.h>
#include <rk3.h>

#include <mpi_plan.h>
#include <tensor.h>
#include <decomp.h>
#include <transpose.h>



void chanfast4(){
  // validateConfig();
  
  initMpiPlan();
  initDecompPlan();
  initFftPlan();
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
    for (int step = 0; step < 3; step++){
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
  // initIo();
  // initGhostzPlan();
  // initTransposePlan();
  // initFftwPlan();
  // initDecompPlan();
  // initMpiPlan();

}