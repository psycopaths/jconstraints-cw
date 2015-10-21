/*
 * Copyright (C) 2015, United States Government, as represented by the 
 * Administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 *
 * The PSYCO: A Predicate-based Symbolic Compositional Reasoning environment 
 * platform is licensed under the Apache License, Version 2.0 (the "License"); you 
 * may not use this file except in compliance with the License. You may obtain a 
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 * Unless required by applicable law or agreed to in writing, software distributed 
 * under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the 
 * specific language governing permissions and limitations under the License.
 */
package gov.nasa.jpf.constraints.solvers.cw;

import java.util.Properties;

import gov.nasa.jpf.constraints.api.ConstraintSolver;
import gov.nasa.jpf.constraints.solvers.ConstraintSolverProvider;
import gov.nasa.jpf.constraints.solvers.cw.CWSolver.CWSolverBuilder;
import gov.nasa.jpf.constraints.solvers.cw.exceptions.CWConfigurationException;

public class CWSolverProvider implements ConstraintSolverProvider {

  @Override
  public String[] getNames() {
    return new String[]{"cwsolver"};
  }

  @Override
  public ConstraintSolver createSolver(Properties config) {
    CWSolverBuilder solverBuilder = new CWSolverBuilder();
    solverBuilder.setConfig(config);
    try {
      if(config.containsKey(CWConfig.LINEAR_SOLVER.getPropStr())) {
        String linearSolverStr = config.getProperty(CWConfig.LINEAR_SOLVER.getPropStr());
        solverBuilder.setLinearSolver(linearSolverStr);
      }
      if(config.containsKey(CWConfig.ITERATIONS_PER_CONSTRAINT.getPropStr())) {
        solverBuilder.setIterationsPerConstraint(Integer.parseInt(config.getProperty(CWConfig.ITERATIONS_PER_CONSTRAINT.getPropStr())));
      }
      if(config.containsKey(CWConfig.NEIGHBORS_GENERATED_PER_ITERATION.getPropStr())) {
        solverBuilder.setNeighborsGeneratedPerIteration(Integer.parseInt(config.getProperty(CWConfig.NEIGHBORS_GENERATED_PER_ITERATION.getPropStr())));
      }
      if(config.containsKey(CWConfig.TABU_ITERATIONS_PER_VARIABLE.getPropStr())) {
        solverBuilder.setTabuIterationsPerVariable(Float.parseFloat(config.getProperty(CWConfig.TABU_ITERATIONS_PER_VARIABLE.getPropStr())));
      }
      if(config.containsKey(CWConfig.MIN_TABU_ITERATIONS.getPropStr())) {
        solverBuilder.setMinTabuIterations(Integer.parseInt(config.getProperty(CWConfig.MIN_TABU_ITERATIONS.getPropStr())));
      }
      if(config.containsKey(CWConfig.ENABLE_SEEDING.getPropStr())) {
        solverBuilder.setEnableSeeding(Boolean.parseBoolean(config.getProperty(CWConfig.ENABLE_SEEDING.getPropStr())));
      }
      if(config.containsKey(CWConfig.ENABLE_BISECTION.getPropStr())) {
        solverBuilder.setEnableBisection(Boolean.parseBoolean(config.getProperty(CWConfig.ENABLE_BISECTION.getPropStr())));
      }
      if(config.containsKey(CWConfig.RANDOMIZATION_RADIUS_EXPONENT.getPropStr())) {
        solverBuilder.setRandomizationRadiusExp(Integer.parseInt(config.getProperty(CWConfig.RANDOMIZATION_RADIUS_EXPONENT.getPropStr())));
      }
      if(config.containsKey(CWConfig.SEED.getPropStr())) {
        solverBuilder.setSeed(Long.parseLong(config.getProperty(CWConfig.SEED.getPropStr())));
      }
      
    } catch(Exception e) {
      throw new CWConfigurationException("Invalid configuration of concolic walk solver.", e);
    }
    return solverBuilder.buildSolver();
  }

}
