/*
 * Problem.cpp
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

#include "Problem.h"
#include <iostream>
#include <sstream>
#include <fstream>

#ifdef WITH_OPENGM
#include <opengm/datastructures/marray/marray_hdf5.hxx> // for hdf5 support
#endif

namespace daoopt {

//extern time_t time_start;

void Problem::removeEvidence() {

  assert(m_n!=UNKNOWN);

  // record original no. of variables
  m_nOrg = m_n;

  // Declare aux. variables.
  int idx, i, new_r, new_n;
  val_t k, new_k;
  //map<unsigned int, unsigned int> evidence(m_evidence);
  vector<val_t> new_domains;
  vector<Function*> new_funs;

  // eliminateVar[i]==TRUE iff var. is to be eliminated
  vector<bool> eliminateVar(m_n,false);
  for (map<int,val_t>::iterator it = m_evidence.begin(); it != m_evidence.end(); ++it) {
    eliminateVar[it->first] = true;
  }

  // Identify and tag unary domain variables
  for (i = 0; i < m_n; ++i) {
    if (m_domains.at(i) == 1) { // regard unary domains as evidence
      m_evidence.insert(make_pair(i, 0));
      ++m_e;
      eliminateVar.at(i) = true;
    }
  }

  // Identify variables not covered by any function
  vector<bool> covered(m_n, false);
  BOOST_FOREACH(Function * f, m_functions) {
    BOOST_FOREACH(int i, f->getScopeVec()) {
      covered.at(i) = true;
    }
  }
  for (size_t i=0; i<m_n; ++i) {
    if (!covered.at(i)) eliminateVar.at(i) = true;
  }

  // Project functions to account for evidence
  m_globalConstant = ELEM_ONE;
  new_r = 0; // max. arity
  vector<Function*>::iterator fi = m_functions.begin();
  for (; fi != m_functions.end(); ++fi) {
    Function *fn = (*fi);
    // Substitute evidence variables
    Function* new_fn = fn->substitute(m_evidence);
    if (new_fn->isConstant()) { // now empty scope
      m_globalConstant OP_TIMESEQ new_fn->getTable()[0];
      delete new_fn;
    } else {
      new_funs.push_back(new_fn); // record new function
      new_r = max(new_r, (int) new_fn->getScopeVec().size());
    }
    delete fn; // delete old function

  }
  m_functions = new_funs;

  // Create dummy function for global constant. Technically the field value needs
  // to be reset to ELEM_ONE, but we keep it around for informational purposes --
  // it should never get used in actual computations.
  double* table1 = new double[1];
  table1[0] = m_globalConstant;
  Function* constFun = new FunctionBayes(m_nOrg, this, set<int>(), table1, 1);
  m_functions.push_back(constFun);

  /*
  // === shrink more by eliminating certain vars ===

  // remember which vars are in which function scope
  map<int, set<Function*> > mapF; // TODO


  // first, build primal graph
  Graph G(m_n);
  for (fi=m_functions.begin(); fi!=m_functions.end(); ++fi)
    G.addClique((*fi)->getScope());

  bool repeatLoop = true;
  while (repeatLoop) {

    for (int i=0; i<m_n; ++i) {
      if (eliminateVar.at(i) || !G.hasNode(i)) continue; // skip this node
    }

    break; // TODO

  }
  */

  // eliminate tagged variables and reorder remaining ones
  idx = new_n = new_k = 0;
  for (int i=0; i<m_n; ++i) {
    if (!eliminateVar.at(i)) {

      m_old2new.insert(make_pair(i,idx));

      k = m_domains.at(i);
      new_domains.push_back(k);
      new_k = max(new_k,k);

      ++idx;
      ++new_n;
    }
  }

  // update variable information
  m_domains = new_domains;
  m_n = new_n;
  m_k = new_k;
#ifndef NO_ASSIGNMENT
//  m_curSolution.resize(m_n,UNKNOWN);
#endif


  // translate scopes of the new functions
  for (fi=m_functions.begin(); fi!=m_functions.end(); ++fi)
    (*fi)->translateScope(m_old2new);

  /*
  cout << "Remapped variables:";
  typedef std::map<int,int>::value_type mtype;
  BOOST_FOREACH( mtype t, m_old2new )
    cout << ' ' << t.first << "->" << t.second;
  cout << endl;
  */

  // update function information
  m_c = m_functions.size();
}


bool Problem::parseOrdering(const string& file, vector<int>& elim) const {

  assert(m_n!=UNKNOWN);

  ifstream inTemp(file.c_str());
  inTemp.close();

  if (inTemp.fail()) { // file not existent yet
    return false;
  }

  igzstream in(file.c_str());

  // ignore first line if there's a pound '#' sign (comment)
  if (in.peek() == '#') {
    in.ignore(8192,'\n');
  }

  int nIn;
  in >> nIn; // length of ordering
  if (nIn != m_n && nIn != m_nOrg) {
    cerr << "Problem reading ordering, number of variables doesn't match" << endl;
    in.close(); exit(1);
  }

  // read into buffer first
  list<int> buffer;
  int x=UNKNOWN;
  while(nIn-- && !in.eof()) {
    in >> x;
    buffer.push_back(x);
  }

  bool fullOrdering = false; // default
  if (buffer.size() == (size_t) m_nOrg) {
    fullOrdering = true;
  }

  int n=0;
  vector<bool> check(m_n, false);

  for (list<int>::iterator it=buffer.begin(); it!=buffer.end(); ++it) {
    if (!fullOrdering) {

      if (*it < 0 || *it >= m_n) {
        cerr << "Problem reading ordering, variable index " << *it << " out of range" << endl;
        in.close(); exit(1);
      }

      if (check[*it]) {
        cerr << "Problem reading ordering, variable " << *it << " appears more than once." << endl;
        in.close(); exit(1);
      } else check[*it] = true;

      elim.push_back(*it); ++n;

    } else { // full order, needs filtering

      if (*it < 0 || *it >= m_nOrg) {
        cerr << "Problem reading ordering, variable index " << *it << " out of range" << endl;
        in.close(); exit(1);
      }

      map<int,int>::const_iterator it2 = m_old2new.find(*it);
      if (it2 != m_old2new.end()) {
        x = it2->second;
        if (check[x]) {
          cerr << "Problem reading ordering, variable " << *it << " appears more than once." << endl;
          in.close(); exit(1);
        } else check[x] = true;

        elim.push_back(x); ++n;
      } else { /* evidence */ }

    }

  }

  if (n!=m_n) {
    cerr << "Problem reading ordering, number of variables doesn't match." << endl;
    in.close();
    exit(1);
  }

  in.close();
  return true;
}


void Problem::saveOrdering(const string& file, const vector<int>& elim) const {
  assert( (int) elim.size() == m_n);

  if (file.substr(file.size()-3,3) == ".gz") { // write gzipped

    ogzstream out(file.c_str());
    if ( ! out ) {
      cerr << "Error writing ordering to file " << file << endl;
      exit(1);
    }
    out << "# daoopt ordering for " << m_name << endl << elim.size();
    for (vector<int>::const_iterator it=elim.begin(); it!=elim.end(); ++it)
      out << ' ' << *it;
    out << endl;
    out.close();

  } else { // write straight text

    ofstream out(file.c_str());
    if ( ! out ) {
      cerr << "Error writing ordering to file " << file << endl;
      exit(1);
    }
    out << "# daoopt ordering for " << m_name << endl << elim.size();
    for (vector<int>::const_iterator it=elim.begin(); it!=elim.end(); ++it)
      out << ' ' << *it;
    out << endl;
    out.close();

  }

}



bool Problem::parseUAI(const string& prob, const string& evid) {
  {
    ifstream inTemp(prob.c_str());
    inTemp.close();

    if (inTemp.fail()) { // file not existent yet
      cerr << "Error reading problem file " << prob << ", aborting." << endl;
      return false;
    }
  }
  if (!evid.empty()) {
    ifstream inTemp(evid.c_str());
    inTemp.close();

    if (inTemp.fail()) { // file not existent yet
      cerr << "Error reading evidence file " << evid << ", aborting." << endl;
      return false;
    }
  }

  igzstream in(prob.c_str());

  // Extract the filename without extension.
  string fname = prob;
  size_t len, start, pos1, pos2;
//  static const basic_string <char>::size_type npos = -1;
#if defined(WINDOWS)
  pos1 = fname.find_last_of("\\");
#elif defined(LINUX)
  pos1 = fname.find_last_of("/");
#endif
  pos2 = fname.find_last_of(".");
  if (pos1 == string::npos) { len = pos2; start = 0; }
  else { len = (pos2-pos1-1); start = pos1+1; }
  m_name = fname.substr(start, len);

  cout << "Reading problem " << m_name << " ..." << endl;

  vector<int> arity;
  vector<vector<int> > scopes;
  string s;
  int x,y;
  val_t xs;
  unsigned int z;

  in >> s; // Problem type
  if (s == "BAYES") {
    m_task = TASK_MAX;
    m_prob = PROB_MULT;
  } else if (s == "MARKOV") {
    m_task = TASK_MAX;
    m_prob = PROB_MULT;
  } else {
    cerr << "Unsupported problem type \"" << s << "\", aborting." << endl;
    in.close(); return false;
  }

  in >> x; // No. of variables
  m_n = x;
  m_domains.resize(m_n,UNKNOWN);
#ifndef NO_ASSIGNMENT
//  m_curSolution.resize(m_n,UNKNOWN);
#endif
  m_k = -1;
  for (int i=0; i<m_n; ++i) { // Domain sizes
    in >> x; // read into int first
    if (x > numeric_limits<val_t>::max()) {
      cerr << "Domain size " << x << " out of range for internal representation.\n"
           << "(Recompile with different type for variable values.)" << endl;
      in.close(); return false;
    }
    xs = (val_t)x;
    m_domains[i] = xs;
    m_k = max(m_k,xs);
  }

  in >> x; // No. of functions
  m_c = x;
  scopes.reserve(m_c);

  // Scope information for functions
  m_r = -1;
  for (int i = 0; i < m_c; ++i)
  {
    vector<int> scope;
    in >> x; // arity

    m_r = max(m_r, x);
    for (int j=0; j<x; ++j) {
      in >> y; // the actual variables in the scope
      if(y>=m_n) {
        cerr << "Variable index " << y << " out of range." << endl;
        in.close(); return false;
      }
      scope.push_back(y); // preserve order from file
    }
    scopes.push_back(scope);
  }

  // Read functions
  for (int i = 0; i < m_c; ++i)
  {
    in >> z; // No. of entries
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
    for (size_t j=0; j<tab_size; ++j) {
      in >> temp[j];
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
  in.close();

  // Read evidence?
  if (evid.empty()) {
    m_e = 0;
    return true; // No evidence, return
  }

  cout << "Reading evidence..." << endl;

  in.open(evid.c_str());

  in >> x;
  m_e = x; // Number of evidence

  for (int i=0; i<m_e; ++i) {
    in >> x; // Variable index
    in >> y; // Variable value
    xs = (val_t) y;
    if (xs >= m_domains[x]) {
      cout << "Variable " << x << " has domain size " << (int) m_domains[x]
           << ", evidence value " << y << " out of range." << endl;
      in.close(); return false;
    }
    m_evidence.insert(make_pair(x,xs));
  }

  in.close();
  return true;
}


void Problem::outputAndSaveSolution(const string& file, const SearchStats* nodestats,
    const vector<count_t>& nodeProf, const vector<count_t>& leafProf, bool toScreen) const {

  bool writeFile = false;
  if (! file.empty())
    writeFile = true;

#ifdef WITH_OPENGM
  // If WITH_OPENGM is defined and file ends with ".h5" the hdf5 output will be used insted of daoopt's own file format.
  // Only the optimal value and its corresponding argument will be stored.

  // chech file ending
  bool useHDF5 = false;
  if(file.size() > 3) {
	  if(file.substr(file.size() - 3) == ".h5") {
	    useHDF5 = true;
	  }
  }

  if(useHDF5) {
	  // disable daoopt's own file output
    writeFile = false;

    // create file
    hid_t handle = marray::hdf5::createFile(file);

    // store value
    vector<double> value(1, m_curCost); // wrapper to use opengm build in hdf5 support
    marray::hdf5::save(handle, "value", value);

    // store argument

    // TODO marray::hdf5::save does not support type signed char this will be fixed in opengm version 2.1.3
    // marray::hdf5::save(handle, "states", m_curSolution);

    // typecast of curSolution to support opengm versions below 2.1.3
    // can be removed when the new opengm version is released.
    vector<size_t> states(m_curSolution.begin(), m_curSolution.end());
    marray::hdf5::save(handle, "states", states);
  }
#endif

  ogzstream out;
  if (writeFile) {
    out.open(file.c_str(), ios::out | ios::trunc | ios::binary);
    if (!out) {
      cerr << "Error writing optimal solution to file " << file << endl;
      writeFile = false;
      out.close();
    }
  }

  oss screen;
  screen << "s " << SCALE_LOG(m_curCost);
#ifndef NO_ASSIGNMENT
  int32_t assigSize = UNKNOWN;
  if (m_subprobOnly)
    assigSize = (int32_t) m_curSolution.size();  // no dummy variable included
  else
    assigSize = (int32_t) m_nOrg;
  screen << ' ' << assigSize;
#endif

  if (writeFile) {
    BINWRITE(out, m_curCost); // mpe solution cost
    count_t countOR = 0, countAND = 0;
    if (nodestats) {
      countOR = nodestats->numExpOR;
      countAND = nodestats->numExpAND;
    }
    BINWRITE(out, countOR);
    BINWRITE(out, countAND);
  }

#ifndef NO_ASSIGNMENT
  if (writeFile) {
    BINWRITE(out,assigSize); // no. of variables in opt. assignment
  }

  // generate full assignment (incl. evidence) if necessary and output
  vector<val_t> outputAssig;
  assignmentForOutput(outputAssig);
  BOOST_FOREACH( int32_t v, outputAssig ) {
    screen << ' ' << v;
    if (writeFile) BINWRITE(out, v);
  }
#endif

  // output node profiles in case of subproblem processing
  if (m_subprobOnly) {
    int32_t size = (int32_t) leafProf.size();
    BINWRITE(out, size);
    // leaf nodes first
    for (vector<count_t>::const_iterator it=leafProf.begin(); it!=leafProf.end(); ++it) {
      BINWRITE(out,*it);
    }
    // now full node profile (has same array size)
    for (vector<count_t>::const_iterator it=nodeProf.begin(); it!=nodeProf.end(); ++it) {
      BINWRITE(out,*it);
    }
  }

  screen << endl;
  if (toScreen)
    cout << screen.str();
  if (writeFile)
    out.close();
}


#ifndef NO_ASSIGNMENT
void Problem::assignmentForOutput(vector<val_t>& assg) const {
  assignmentForOutput(m_curSolution, assg);
}

void Problem::assignmentForOutput(const vector<val_t>& inAssg, vector<val_t>& outAssg) const {
  if (m_subprobOnly || inAssg.empty()) {
    outAssg = inAssg;
    // update: no need to remove dummy anymore
  } else {
    outAssg.resize(m_nOrg, UNKNOWN);
    for (int i=0; i<m_nOrg; ++i) {
      map<int,int>::const_iterator itRen = m_old2new.find(i);
      if (itRen != m_old2new.end()) {  // var part of solution
        outAssg.at(i) = inAssg.at(itRen->second);
      } else {
        map<int,val_t>::const_iterator itEvid = m_evidence.find(i);
        if (itEvid != m_evidence.end())  // var part of evidence
          outAssg.at(i) = itEvid->second;
        else  // var had unary domain
          outAssg.at(i) = 0;
      }
    }
  }
}
#endif


void Problem::updateSolution(double cost,
#ifndef NO_ASSIGNMENT
    const vector<val_t>& sol,
#endif
    const SearchStats* nodestats,
    bool output) {

  if (ISNAN(cost))
    return;

  double costCheck = ELEM_ZERO;
#ifndef NO_ASSIGNMENT
  // check for complete assignment first
  for (size_t i = 0; i < sol.size(); ++i) {
    if (sol[i] == NONE) {
      oss ss; ss << "Warning: skipping incomplete solution, reported " << cost;
      DIAG(ss << " " << sol.size() << " " << sol;)
      ss << endl; myprint(ss.str());
      return;
    }
    if (!m_subprobOnly && sol[i] >= m_domains[i]) {
      oss ss; ss << "Warning: value " << (int)sol[i] << " outside of variable " << i
                 << " domain " << (int)m_domains[i];
      ss << endl; myprint(ss.str());
      return;
    }
  }

  // use Kahan summation to compute exact solution cost
  // TODO (might not work in non-log scale?)
  if (cost != ELEM_ZERO && !m_subprobOnly) {
    costCheck = ELEM_ONE; double comp = ELEM_ONE;  // used across loop iterations
    double y, z;  // reset for each loop iteration
    BOOST_FOREACH( Function* f, m_functions ) {
      z = f->getValue(sol);
//      if (z == ELEM_ZERO) {
//        oss ss; ss << "Warning: skipping zero-cost solution. Reported cost: " << cost;
//        DIAG(ss << " " << sol.size() << " " << sol;)
//        ss << endl; myprint(ss.str());
//        return;
//      }
      y = z OP_DIVIDE comp;
      z = costCheck OP_TIMES y;
      comp = (z OP_DIVIDE costCheck) OP_DIVIDE y;
      costCheck = z;
    }
    if (1e-3 < abs(cost-costCheck)) {
      oss ss; ss << "Warning: solution cost " << costCheck << " differs significantly"
		 << ", reported " << cost << endl;
      myprint(ss.str());
    }
  } else
#endif
  costCheck = cost;

  if (ISNAN(costCheck) || (!ISNAN(m_curCost) && costCheck <= m_curCost)) { // TODO costCheck =?= ELEM_ZERO )
    oss ss; ss << "Warning: Discarding solution with cost " << costCheck << ", reported: " << cost;
#ifndef NO_ASSIGNMENT
    DIAG(ss << " " << sol.size() << " " << sol;)
#endif
    ss << endl; myprint(ss.str());
    return;
  }
  m_curCost = costCheck;
  if (costCheck == ELEM_ZERO) output = false;
  ostringstream ss;
  if (output) {
    ss << "u ";
    if (nodestats)
      ss << nodestats->numExpOR << ' ' <<  nodestats->numExpAND << ' ';
    else
      ss << "0 0 ";
    ss << SCALE_LOG(costCheck) ;
  }

#ifndef NO_ASSIGNMENT
  // save only the reduced solution
  // NOTE: sol.size() < m_nOrg in conditioned subproblem case
  m_curSolution = sol;
  // output the complete assignment (incl. evidence)
  if (output) {
    vector<val_t> outputAssg;
    assignmentForOutput(outputAssg);
    ss << ' ' << outputAssg.size();
    BOOST_FOREACH( int v, outputAssg ) {
      ss << ' ' << v;
    }
  }
#endif

  if (output) {
    ss << endl;
    myprint(ss.str());
  }
}


void Problem::resetSolution() {
  m_curCost = ELEM_NAN;
#ifndef NO_ASSIGNMENT
  m_curSolution.clear();
#endif
}


void Problem::writeUAI(const string& prob) const {
  assert (prob.size());

  ogzstream out;
  out.open(prob.c_str(), ios::out | ios::trunc);

  if (!out) {
    cerr << "Error writing reduced network to file " << prob << endl;
    exit(1);
  }

  out << "MARKOV" << endl; // TODO hard-coding is not optimal

  // variable info
  out << (m_n) << endl;
  for (vector<val_t>::const_iterator it=m_domains.begin(); it!=m_domains.end(); ++it)
    out << ' ' << ((int) *it) ;

  // function information
  out << endl << m_functions.size() << endl;
  for (vector<Function*>::const_iterator it=m_functions.begin(); it!=m_functions.end(); ++it) {
    const vector<int>& scope = (*it)->getScopeVec();
    out << scope.size() << '\t'; // scope size
    for (vector<int>::const_iterator itS=scope.begin(); itS!=scope.end(); ++itS)
      out << *itS << ' '; // variables in scope
    out << endl;
  }
  out << endl;

  // write the function tables
  for (vector<Function*>::const_iterator it=m_functions.begin(); it!=m_functions.end(); ++it) {
    double * T = (*it)->getTable();
    out << (*it)->getTableSize() << endl; // table size
    for (size_t i=0; i<(*it)->getTableSize(); ++i)
      out << ' ' << SCALE_NORM( T[i] ); // table entries
    out << endl;
  }

  // done
  out << endl;
  out.close();

}


void Problem::addDummy() {
  m_n += 1;
  m_hasDummy = true;
  m_domains.push_back(1); // unary domain
}


void Problem::replaceFunctions(const vector<Function*>& newFunctions) {
  // delete current functions
  for (vector<Function*>::iterator it = m_functions.begin(); it!= m_functions.end(); ++it) {
    if (*it) delete (*it);
  }
  // store new functions
  m_functions = newFunctions;
  m_c = m_functions.size();
  // update function scopes???
}


#ifndef NO_ASSIGNMENT
bool Problem::isEliminated(int i) const {
  map<int,int>::const_iterator itRen = m_old2new.find(i);
  return itRen == m_old2new.end();
}
#endif

}  // namespace daoopt
