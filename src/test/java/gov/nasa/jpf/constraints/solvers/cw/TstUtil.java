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
import gov.nasa.jpf.constraints.api.ConstraintSolver.Result;
import gov.nasa.jpf.constraints.api.Expression;
import gov.nasa.jpf.constraints.api.Valuation;
import gov.nasa.jpf.constraints.api.ValuationEntry;
import gov.nasa.jpf.constraints.solvers.ConstraintSolverFactory;
import gov.nasa.jpf.constraints.solvers.cw.CWSolver.CWSolverBuilder;

public class TstUtil {
  
	public static CWSolver createCWSolver(Properties conf) {
		conf.setProperty("symbolic.dp", "cwsolver");
		ConstraintSolverFactory factory = new ConstraintSolverFactory(conf);
		ConstraintSolver solver = factory.createSolver();
		return (CWSolver) solver;
	}
	
	public static CWSolver createCWSolver() {
	  CWSolverBuilder builder = new CWSolverBuilder();
	  builder.setIterationsPerConstraint(300);
    return builder.buildSolver();
  }
	
	public static Valuation runTest(Expression<Boolean> expr, Result expectedRes, boolean printExpr) {
	  return runTest(createCWSolver(), expr, expectedRes, printExpr);
	}
	
	 public static Valuation runTest(Expression<Boolean> expr, Result expectedRes) {
	    return runTest(createCWSolver(), expr, expectedRes, true);
	  }
	
	public static Valuation runTest(ConstraintSolver solver, Expression<Boolean> expr, Result expectedRes, boolean printExpr) {
		if(printExpr)
		  print("Expr: " + expr.toString());
		try {
			Valuation val = new Valuation();
			long start = System.currentTimeMillis();
	        Result res = solver.solve(expr, val);
	        long solverTime = (System.currentTimeMillis() - start);
	        print("Solver time: " + solverTime + "ms");
	        print("Expected " + expectedRes + " got " + res);
	        if(res == Result.SAT) {
	          print("-------Valuation-------");
		        for(ValuationEntry<?> exp : val)
		          print(exp.getVariable() + "=" + exp.getValue());
		        print("-----------------------");
	        }
	        junit.framework.Assert.assertEquals(res, expectedRes);
	        return val;
		} catch(Exception e) {
			throw e;
		} finally {
		  print("======================================================================");
		}
	}
	
	public static void print(String str) {
	  System.out.println("[Junit] " + str);
	}

}
