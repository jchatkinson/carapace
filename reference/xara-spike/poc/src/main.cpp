// Xara wasm POC: hand-wired analysis chain (bypassing BasicAnalysisBuilder /
// ModelRegistry, which eagerly pull in XaraClassBroker's ~230-class catalog).
// Wiring sequence copied from SRC/runtime/runtime/BasicAnalysisBuilder.cpp
// (initialize/number/domainChanged/analyzeStatic), since that's the real,
// tested composition logic -- just invoked directly instead of through the
// class that also drags in the full element/material broker.
//
// Model: single horizontal 2-node truss, node 1 fixed, node 2 free & loaded
// axially. Known answer: displacement = P*L / (A*E).

#include <Domain.h>
#include <Node.h>
#include <Truss.h>
#include <ElasticMaterial.h>
#include <SP_Constraint.h>
#include <StaticPattern.h>
#include <NodalLoad.h>
#include <LinearSeries.h>

#include <AnalysisModel.h>
#include <PlainHandler.h>
#include <RCM.h>
#include <DOF_Numberer.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#include <LoadControl.h>
#include <Linear.h>
#include <CTestNormUnbalance.h>

#include <cstdio>
#include <cmath>

int main() {
  const double L = 100.0;
  const double A = 2.0;
  const double E = 30000.0;
  const double P = 50.0;

  Domain domain;

  Node *n1 = new Node(1, 2, 0.0, 0.0);
  Node *n2 = new Node(2, 2, L,   0.0);
  domain.addNode(n1);
  domain.addNode(n2);

  domain.addSP_Constraint(new SP_Constraint(1, 0, 0.0, true));
  domain.addSP_Constraint(new SP_Constraint(1, 1, 0.0, true));
  domain.addSP_Constraint(new SP_Constraint(2, 1, 0.0, true));

  UniaxialMaterial *mat = new ElasticMaterial(1, E);
  Truss *truss = new Truss(1, 2, 1, 2, *mat, A);
  domain.addElement(truss);

  LinearSeries *series = new LinearSeries();
  StaticPattern *pattern = new StaticPattern(1, 1.0, series);
  domain.addLoadPattern(pattern);

  static double loadValues[2] = {P, 0.0};
  Vector loadVec(loadValues, 2);
  NodalLoad *load = new NodalLoad(1, 2, loadVec);
  domain.addNodalLoad(load, 1);

  AnalysisModel      theModel(domain);
  PlainHandler       theHandler;
  DOF_Numberer       theNumberer(*(new RCM(false)));
  BandGenLinSOE       theSOE(*(new BandGenLinLapackSolver()));
  LoadControl        theIntegrator(1.0, 1, 1.0, 1.0);
  Linear             theAlgorithm;
  CTestNormUnbalance theTest(1.0e-8, 25, 0);

  theModel.setLinks(domain, theHandler);
  theHandler.setLinks(theModel);
  theIntegrator.setLinks(theModel, theSOE, &theTest);
  theAlgorithm.setLinks(theIntegrator, theSOE, &theTest);

  // ---- number() ----
  theModel.clearAll();
  if (theHandler.handle() < 0) { printf("FAIL: handle()\n"); return 1; }
  theNumberer.setLinks(theModel);
  if (theNumberer.numberDOF() < 0) { printf("FAIL: numberDOF()\n"); return 1; }
  if (theHandler.doneNumberingDOF() < 0) { printf("FAIL: doneNumberingDOF()\n"); return 1; }

  // ---- domainChanged() ----
  Graph &theGraph = theModel.getDOFGraph();
  if (theSOE.setSize(theGraph) < 0) { printf("FAIL: SOE.setSize()\n"); return 1; }
  theModel.clearDOFGraph();
  if (theIntegrator.domainChanged() < 0) { printf("FAIL: integrator.domainChanged()\n"); return 1; }

  // ---- initialize() ----
  if (theIntegrator.initialize() < 0) { printf("FAIL: integrator.initialize()\n"); return 1; }
  theIntegrator.commit();
  domain.initialize();

  // ---- analyzeStatic(1 step) ----
  if (theModel.analysisStep(0.0) < 0) { printf("FAIL: analysisStep()\n"); return 1; }
  if (theIntegrator.newStep() < 0) { printf("FAIL: newStep()\n"); return 1; }
  if (theAlgorithm.solveCurrentStep() < 0) { printf("FAIL: solveCurrentStep()\n"); return 1; }
  if (theIntegrator.commit() < 0) { printf("FAIL: commit()\n"); return 1; }

  double dx = n2->getDisp()(0);
  double expected = P * L / (A * E);

  printf("computed dx = %.10f, expected = %.10f\n", dx, expected);
  bool pass = std::fabs(dx - expected) < 1e-6;
  printf(pass ? "PASS\n" : "FAIL\n");
  return pass ? 0 : 1;
}
