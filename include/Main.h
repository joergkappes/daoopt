/*
 * Main.h
 *
 *  Copyright (C) 2008-2012 Lars Otten
 *  This file is part of DAOOPT.
 *
 *  DAOOPT is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DAOOPT is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DAOOPT.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  Created on: Oct 18, 2011
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef MAIN_H_
#define MAIN_H_

#include "_base.h"

#include "Problem.h"
#include "Function.h"
#include "Graph.h"
#include "Pseudotree.h"
#include "ProgramOptions.h"
#include "MiniBucketElim.h"
#ifdef ENABLE_SLS
#include "SLSWrapper.h"
#endif

#ifdef PARALLEL_DYNAMIC
  #include "BranchAndBoundMaster.h"
  #include "BoundPropagatorMaster.h"
  #include "SubproblemHandler.h"
  #include "SigHandler.h"
#else
  #ifdef PARALLEL_STATIC
    #include "ParallelManager.h"
    #include "BranchAndBoundSampler.h"
  #endif
  #include "BranchAndBound.h"
  #include "BranchAndBoundRotate.h"
  #include "BoundPropagator.h"
#endif

#include "BestFirst.h"
#include "LimitedDiscrepancy.h"

namespace daoopt {

class Main {
protected:
  bool m_solved;
  scoped_ptr<ProgramOptions> m_options;
  scoped_ptr<Problem> m_problem;
  scoped_ptr<Pseudotree> m_pseudotree;
  scoped_ptr<Heuristic> m_heuristic;
#ifdef ENABLE_SLS
  scoped_ptr<SLSWrapper> m_slsWrapper;
#endif

#if defined PARALLEL_DYNAMIC
  scoped_ptr<BranchAndBoundMaster> m_search;
  scoped_ptr<SearchSpaceMaster> m_space;
#else
#ifdef PARALLEL_STATIC
  scoped_ptr<ParallelManager> m_search;
#else
  scoped_ptr<Search> m_search;
#endif
  scoped_ptr<SearchSpace> m_space;
#endif

protected:
  bool runSearchDynamic();
  bool runSearchStatic();
  bool runSearchWorker();

public:
  bool start() const;
  bool parseOptions(int argc, char** argv);
  bool outputInfo() const;
  bool loadProblem();
  bool findOrLoadOrdering();
  bool runSLS();
  bool initDataStructs();
  bool compileHeuristic();
  bool runLDS();
  bool finishPreproc();
  bool runSearch();
  bool outputStats() const;

  bool isSolved() const { return m_solved; }

  Main();

  // getter / setter
  void setOptions(ProgramOptions* options); // takes ownership of input
  void setProblem(Problem* problem); // takes ownership of input
  const Problem& getProblem() const;

};

/* Inline implementations */

inline Main::Main() : m_solved(false) {
  /* nothing here */
}

inline bool Main::runSearch() {
  if (m_options->nosearch)
    return true;
  if (m_solved)
    return true;
  cout << "--- Starting search ---" << endl;
#if defined PARALLEL_DYNAMIC
  return runSearchDynamic();
#elif defined PARALLEL_STATIC
  return runSearchStatic();
#else
  return runSearchWorker();
#endif
}

}  // namespace daoopt

#endif /* MAIN_H_ */
