/*
 * Problem.h
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
 *  Created on: Oct 10, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "Function.h"
#include "Graph.h"
#include "MiniBucket.h"
#include "SearchSpace.h"
#include "_base.h"
#include "gzstream.h"

#ifdef WITH_OPENGM
#include <opengm/utilities/indexing.hxx>
#endif
namespace daoopt {

/* holds a problem instance with variable domains and function tables */
class Problem {

protected:

  bool m_is_copy;        // true iff object instance is copy

  bool m_subprobOnly;    // Solving only a conditioned subproblem
  bool m_hasDummy;       // is last variable a dummy variable?

  int m_prob;            // Problem class (multiplication or summation of costs)
  int m_task;            // Type of problem (Minim. or maxim. task)


  int m_n;               // No. of variables
  int m_nOrg;            // No. of variables (before evidence was removed)
  val_t m_k;             // Max. domain size
  int m_e;               // No. of evidence

  int m_c;               // No. of functions
  int m_r;               // Max. function arity

  double m_globalConstant;        // Global constant modifier for objective function

  double m_curCost;               // Cost of current solution

  string m_name;                  // Problem name

  vector<val_t> m_domains;        // Domain sizes of variables

  vector<Function*> m_functions;  // List of functions

  map<int,val_t> m_evidence;      // List of evidence as <index,value>

  map<int,int> m_old2new;         // Translation of variable names after removing evidence

  vector<val_t> m_curSolution;       // Current best solution

public:
  void setCopy() { m_is_copy = true; }

  val_t getDomainSize(int i) const;
  double globalConstInfo() const;

  int getN() const { return m_n; }
  int getNOrg() const { return m_nOrg; }
  val_t getK() const { return m_k; }
  int getE() const { return m_e; }
  int getC() const { return m_c; }
  int getR() const { return m_r; }

  void setSubprobOnly() { m_subprobOnly = true; }
  const string& getName() const { return m_name; }

  const vector<Function*>& getFunctions() const { return m_functions; }
  const vector<val_t>& getDomains() const { return m_domains; }

  // replaces the current set of functions with an equivalent one
  // (pseudo tree compatibility is implicitly assumed)
  void replaceFunctions(const vector<Function*>& newFunctions);

  bool hasDummy() const { return m_hasDummy; }

public:

  /* parses a UAI format input file */
  bool parseUAI(const string& prob, const string& evid);

#ifdef WITH_OPENGM
  /* convert OpenGM model */
  template <class GM>
  bool convertOPENGM(const GM& gm);
#endif

  /* writes the current problem to a UAI file */
  void writeUAI(const string& prob) const;

  /* parses an ordering from file 'file' and stores it in 'elim' */
  bool parseOrdering(const string& file, vector<int>& elim) const;
  /* stores ordering from 'elim' in file 'file' */
  void saveOrdering(const string& file, const vector<int>& elim) const;

  /* removes evidence and unary-domain variables */
  void removeEvidence();

  /* retrieve the current optimal solution */
  double getSolutionCost() const { return m_curCost; }

#ifndef NO_ASSIGNMENT
  /* retrieve the current optimal assignment */
  const vector<val_t>& getSolutionAssg() const { return m_curSolution; }
  /* compute the assignment for output (might add evidence back in) */
  void assignmentForOutput(vector<val_t>&) const;
  void assignmentForOutput(const vector<val_t>& in, vector<val_t>& out) const;
#endif

  /* report an updated solution */
  void updateSolution(double cost,
#ifndef NO_ASSIGNMENT
      const vector<val_t>& sol,
#endif
      const SearchStats* nodestats = NULL,
      bool output = true);

  /* resets current optimal solution cost and assignment */
  void resetSolution();

  /* outputs the solution to the screen and, if file!="", writes it to file
   * (for subproblem solving, only relevant variables will be output)
   *  - cost is the MPE tuple value
   *  - sol is the optimal solution tuple
   *  - noNodes is the number of OR/AND nodes
   *  - nodeProf and leafProf are the full and leaf node profiles
   *  - if toScreen==true, will skip the console output (file only)
   */
  void outputAndSaveSolution(const string& file, const SearchStats* nodestats,
                             const vector<count_t>& nodeProf, const vector<count_t>& leafProf,
                             bool toScreen = true) const;

#ifndef NO_ASSIGNMENT
  /* returns true iff the index variable from the full set has been eliminated
   * as evidence or unary */
  bool isEliminated(int i) const;
#endif

  /* adds the dummy variable to connect disconnected pseudo tree components */
  void addDummy();

public:
  Problem();
  virtual ~Problem();
};


/* Inline definitions */

inline val_t Problem::getDomainSize(int i) const {
  assert (i<m_n);
  return m_domains[i];
}

inline double Problem::globalConstInfo() const {
  return m_globalConstant;
}

inline Problem::Problem() :
    m_is_copy(false),
    m_subprobOnly(false),
    m_hasDummy(false),
    m_prob(UNKNOWN),
    m_task(UNKNOWN),
    m_n(UNKNOWN),
    m_nOrg(UNKNOWN),
    m_k(UNKNOWN),
    m_e(UNKNOWN),
    m_c(UNKNOWN),
    m_r(UNKNOWN),
    m_globalConstant(ELEM_NAN),
    m_curCost(ELEM_NAN)
{ /* empty*/ }

inline Problem::~Problem() {
  // delete functions
  if (!m_is_copy)
    for (vector<Function*>::iterator it = m_functions.begin(); it!= m_functions.end(); ++it)
      if (*it) delete (*it);
}

#ifdef WITH_OPENGM
/* convert OpenGM model */
template <class GM>
inline bool Problem::convertOPENGM(const GM& gm) {
  typedef typename GM::FactorType::ShapeIteratorType FactorShapeIteratorType;

  vector<int> arity;
  vector<vector<int> > scopes;
  string s;
  int x,y;
  val_t xs;
  unsigned int z;

  m_task = TASK_MAX;
  m_prob = PROB_MULT;

  m_n = gm.numberOfVariables();
  m_domains.resize(m_n,UNKNOWN);

  m_k = -1;
  for (int i=0; i<m_n; ++i) { // Domain sizes
	x = gm.numberOfLabels(i);
	if (x > numeric_limits<val_t>::max()) {
	 cerr << "Domain size " << x << " out of range for internal representation.\n"
		  << "(Recompile with different type for variable values.)" << endl;
	}
	xs = (val_t)x;
	m_domains[i] = xs;
	m_k = max(m_k,xs);
  }
  x = gm.numberOfFactors();
  m_c = x;
  scopes.reserve(m_c);

  // Scope information for functions
  m_r = -1;
  for (int i = 0; i < m_c; ++i)
  {
	vector<int> scope;
	x = gm[i].numberOfVariables(); // arity

	m_r = max(m_r, x);
	for (int j=0; j<x; ++j) {
	  y=gm[i].variableIndex(j); // the actual variables in the scope
	  if(y>=m_n) {
		cerr << "Variable index " << y << " out of range." << endl;
	  }
		scope.push_back(y); // preserve order from file
	  }
	  scopes.push_back(scope);
  }

  // Read functions
  for (int i = 0; i < m_c; ++i)
  {
	z=gm[i].size(); // No. of entries
	size_t tab_size = 1;

	for (vector<int>::iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it) {
	  tab_size *= m_domains[*it];
	}

	assert(tab_size==z); // product of domain sizes matches no. of entries

	// create set version of the scope (ordered)
	set<int> scopeSet(scopes[i].begin(), scopes[i].end());
	z = scopeSet.size();

	// compute reindexing map from specified scope to ordered, internal one
	map<int,int> mapping; int k=0;
	for (vector<int>::const_iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it)
	  mapping[*it] = k++;
	vector<int> reidx(z);
	vector<int>::iterator itr = reidx.begin();
	for (set<int>::iterator it=scopeSet.begin(); itr!=reidx.end(); ++it, ++itr) {
	  *itr = mapping[*it];
	}

	// read the full table into an temp. array (to allow reordering)
	vector<double> temp(tab_size);
	opengm::ShapeWalkerSwitchedOrder<FactorShapeIteratorType> shapeWalker(gm[i].shapeBegin(),gm[i].numberOfVariables());
	//opengm::ShapeWalker<FactorShapeIteratorType> shapeWalker(gm[i].shapeBegin(),gm[i].numberOfVariables());
	for(typename GM::IndexType scalarIndex=0;scalarIndex<tab_size;++scalarIndex,++shapeWalker) {
	  //temp[scalarIndex]=0.1+scalarIndex;
	  temp[scalarIndex]=pow(10.0,-gm[i](shapeWalker.coordinateTuple().begin()));

	  //temp[scalarIndex]=gm[i](shapeWalker.coordinateTuple().begin());
	}

	// get the variable domain sizes
	vector<val_t> limit; limit.reserve(z);
	for (vector<int>::const_iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it)
	  limit.push_back(m_domains[*it]);
	vector<val_t> tuple(z, 0);

	// create the new table (with reordering)
	double* table = new double[tab_size];
	for (size_t j=0; j<tab_size; ) {
	  size_t pos=0, offset=1;
	  // j is the index in the temp. table
	  for (k=z-1; k>=0; --k) { // k goes backwards through the ordered scope
		pos += tuple[reidx[k]] * offset;
		offset *= m_domains[scopes[i][reidx[k]]];
	  }
	  table[pos] = ELEM_ENCODE( temp[j] );
	  increaseTuple(j,tuple,limit);
	}

	Function* f = new FunctionBayes(i,this,scopeSet,table,tab_size);
	m_functions.push_back(f);

  } // All function tables read
  return true;
}
#endif

}  //

#endif /* PROBLEM_H_ */

