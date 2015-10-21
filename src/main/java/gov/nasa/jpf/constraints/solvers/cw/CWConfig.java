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

import gov.nasa.jpf.constraints.api.ConstraintSolver;
import gov.nasa.jpf.constraints.solvers.ConstraintSolverFactory;

public class CWConfig<T> {
  public static final CWConfig<String> LINEAR_SOLVER = new CWConfig<String>("concolic_walk.linear_solver","z3");
  public static final CWConfig<Integer> ITERATIONS_PER_CONSTRAINT = new CWConfig<Integer>("concolic_walk.iterations_per_constraint", 150);
  public static final CWConfig<Integer> NEIGHBORS_GENERATED_PER_ITERATION = new CWConfig<Integer>("concolic_walk.neighbors_generated_per_iteration", 10);
  public static final CWConfig<Float> TABU_ITERATIONS_PER_VARIABLE = new CWConfig<Float>("concolic_walk.tabu_iterations_per_variable", 0.5f);
  public static final CWConfig<Integer> MIN_TABU_ITERATIONS = new CWConfig<Integer>("concolic_walk.min_tabu_iterations", 3);
  public static final CWConfig<Boolean> ENABLE_SEEDING = new CWConfig<Boolean>("concolic_walk.enable_seeding", true);
  public static final CWConfig<Boolean> ENABLE_BISECTION = new CWConfig<Boolean>("concolic_walk.enable_bisection", true);
  public static final CWConfig<Integer> RANDOMIZATION_RADIUS_EXPONENT = new CWConfig<Integer>("concolic_walk.randomization_radius_exponent", 12);
  public static final CWConfig<Long> SEED = new CWConfig<Long>("concolic_walk.seed", 233L);
  
  private final String opt;
  private final T defaultVal;
  
  private CWConfig(String opt, T defaultVal) {
    this.opt = opt;
    this.defaultVal = defaultVal;
  }
  
  @Override
  public String toString() {
    return this.opt + "=" + this.defaultVal;
  }
  
  public String getPropStr() {
    return this.opt;
  }
  
  public T getDefaultValue() {
    return this.defaultVal;
  }
}
