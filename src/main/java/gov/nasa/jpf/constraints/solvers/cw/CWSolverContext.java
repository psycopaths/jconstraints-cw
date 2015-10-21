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

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import gov.nasa.jpf.constraints.api.ConstraintSolver.Result;
import gov.nasa.jpf.constraints.api.Expression;
import gov.nasa.jpf.constraints.api.SolverContext;
import gov.nasa.jpf.constraints.api.Valuation;
import gov.nasa.jpf.constraints.expressions.LogicalOperator;
import gov.nasa.jpf.constraints.expressions.PropositionalCompound;
import gov.nasa.jpf.constraints.util.ExpressionUtil;

/**
 * TODO: the concolic walk algorithm does NOT support incremental solving yet (even though the underlying linear solver, e.g., z3 does) 
 */
public class CWSolverContext extends SolverContext {

  private CWSolver cwSolver;
  private Deque<List<Expression<Boolean>>> constraints = new ArrayDeque<List<Expression<Boolean>>>();
  
  public CWSolverContext(CWSolver solver) {
    this.cwSolver = solver;
    
    //push the initial context
    this.constraints.push(new LinkedList<Expression<Boolean>>());
  }
  
  @Override
  public Result solve(Valuation val) {
    if(this.constraints.isEmpty())
      return Result.UNSAT;//TODO: Is this the correct semantics?
    /*
     * TODO: maybe not the best way of solving the constraints performance-wise.
     * While processing the stack, we might as well conduct the splitting of
     * linear and non-linear constraints. Anyway, this is easy...
     */
    Expression<Boolean> expr = combineDeque(this.constraints);
    return this.cwSolver.solve(expr, val);
  }
  
  
  private Expression<Boolean> combineDeque(Deque<List<Expression<Boolean>>> deq) {
    List<Expression<Boolean>> cList = new ArrayList<>();
    for(List<Expression<Boolean>> l : deq) {
      cList.add(ExpressionUtil.and(l));
    }
    return ExpressionUtil.and(cList);
  }
  
  @Override
  public void push() {
    constraints.push(new LinkedList<Expression<Boolean>>());
  }
  
  @Override
  public void pop(int n) {
    for(int i = 0; i < n; i++)
      constraints.pop();
  }

  @Override
  public void add(List<Expression<Boolean>> expressions) {
    constraints.peek().addAll(expressions);
  }

  @Override
  public void dispose() {
    // Nothing to do here?
  }
}
