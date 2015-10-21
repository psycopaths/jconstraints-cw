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

import org.junit.Test;

import gov.nasa.jpf.constraints.api.Expression;
import gov.nasa.jpf.constraints.api.Variable;
import gov.nasa.jpf.constraints.api.ConstraintSolver.Result;
import gov.nasa.jpf.constraints.expressions.Constant;
import gov.nasa.jpf.constraints.expressions.LogicalOperator;
import gov.nasa.jpf.constraints.expressions.NumericBooleanExpression;
import gov.nasa.jpf.constraints.expressions.NumericComparator;
import gov.nasa.jpf.constraints.expressions.NumericCompound;
import gov.nasa.jpf.constraints.expressions.NumericOperator;
import gov.nasa.jpf.constraints.expressions.PropositionalCompound;
import gov.nasa.jpf.constraints.expressions.functions.FunctionExpression;
import gov.nasa.jpf.constraints.expressions.functions.math.MathFunctions;
import gov.nasa.jpf.constraints.types.BuiltinTypes;

public class SimpleTest {

  
  @Test
  public void simpleNonlinearConstraint() {
    /*
     * ((('x' * 'x') >= 'x') && (true && (2.42 > 's')))
     */
    Expression<Boolean> expr = new PropositionalCompound(new NumericBooleanExpression(new NumericCompound(new Variable<Double>(BuiltinTypes.DOUBLE, "x"), NumericOperator.MUL, new Variable<Double>(BuiltinTypes.DOUBLE, "x")), NumericComparator.GE, new Variable<Double>(BuiltinTypes.DOUBLE, "x")), LogicalOperator.AND, new PropositionalCompound(new Constant<Boolean>(BuiltinTypes.BOOL, true), LogicalOperator.AND, new NumericBooleanExpression(
        new Constant<Double>(BuiltinTypes.DOUBLE, 2.42),
        NumericComparator.GT,
        new Variable<Double>(BuiltinTypes.DOUBLE, "s"))));
    TstUtil.runTest(expr, Result.SAT);
  }
  
  @Test
  public void simpleNonlinearConstraint2() {
    /*
     * ((('x' * 'x') > 'x') && (true && (2.42 > 's')))
     */
    Expression<Boolean> expr = new PropositionalCompound(new NumericBooleanExpression(new NumericCompound(new Variable<Double>(BuiltinTypes.DOUBLE, "x"), NumericOperator.MUL, new Variable<Double>(BuiltinTypes.DOUBLE, "x")), NumericComparator.GT, new Variable<Double>(BuiltinTypes.DOUBLE, "x")), LogicalOperator.AND, new PropositionalCompound(new Constant<Boolean>(BuiltinTypes.BOOL, true), LogicalOperator.AND, new NumericBooleanExpression(
        new Constant<Double>(BuiltinTypes.DOUBLE, 2.42),
        NumericComparator.GT,
        new Variable<Double>(BuiltinTypes.DOUBLE, "s"))));
    TstUtil.runTest(expr, Result.SAT);
  }
  
  @Test
  public void simpleNonlinearConstraint3() {
    /*
     * ((('x' * 'x') > 'x') && (2.0 > 'x'))
     */
    Expression<Boolean> expr = new PropositionalCompound(new NumericBooleanExpression(new NumericCompound(new Variable<Double>(BuiltinTypes.DOUBLE, "x"), NumericOperator.MUL, new Variable<Double>(BuiltinTypes.DOUBLE, "x")), NumericComparator.GT, new Variable<Double>(BuiltinTypes.DOUBLE, "x")), LogicalOperator.AND, new NumericBooleanExpression(
        new Constant<Double>(BuiltinTypes.DOUBLE, 2.0),
        NumericComparator.GT,
        new Variable<Double>(BuiltinTypes.DOUBLE, "x")));
    TstUtil.runTest(expr, Result.SAT);
  }
  
  
  /*@Test
  public void simpleNonlinearConstraint2() {
    //((('x' * 'x') < 's') && (true && (2.42 > 's')))
    Expression<Boolean> expr = new PropositionalCompound(new NumericBooleanExpression(new NumericCompound(new Variable<Double>(BuiltinTypes.DOUBLE, "x"), NumericOperator.MUL, new Variable<Double>(BuiltinTypes.DOUBLE, "x")), NumericComparator.LT, new Variable<Double>(BuiltinTypes.DOUBLE, "s")), LogicalOperator.AND, new PropositionalCompound(new Constant<Boolean>(BuiltinTypes.BOOL, true), LogicalOperator.AND, new NumericBooleanExpression(
        new Constant<Double>(BuiltinTypes.DOUBLE, 2.42),
        NumericComparator.GT,
        new Variable<Double>(BuiltinTypes.DOUBLE, "s"))));
    TstUtil.runTest(expr, Result.SAT);
  }*/
  
  @Test
  public void simpleTrigConstraint() {
    //Identity:
    //('tan'('x') == ('sin'('x') / 'cos'('x')))
    Variable<Double> x = new Variable<Double>(BuiltinTypes.DOUBLE, "x");
    Expression<Boolean> expr = new NumericBooleanExpression(
        new FunctionExpression<>(MathFunctions.TAN, x),
        NumericComparator.EQ, new NumericCompound<Double>(
            new FunctionExpression<>(MathFunctions.SIN, x),
            NumericOperator.DIV, new FunctionExpression<>(
                MathFunctions.COS, x)));
    
    //Consists of ONLY trig methods (nonlinear).
    //TODO: We cannot solve these with concolic-walk only! We could if complemented by coral
    
    TstUtil.runTest(expr, Result.DONT_KNOW);
  }
  
  @Test
  public void simpleTrigConstraint2() {
    //Identity:
    //(('x' >= 1.2) && ('tan'('x') == ('sin'('x') / 'cos'('x'))))
    Variable<Double> x = new Variable<Double>(BuiltinTypes.DOUBLE, "x");
    Expression<Boolean> expr = new PropositionalCompound(new NumericBooleanExpression(x, NumericComparator.GE, new Constant<Double>(BuiltinTypes.DOUBLE, 1.2)), LogicalOperator.AND, new NumericBooleanExpression(
        new FunctionExpression<>(MathFunctions.TAN, x),
        NumericComparator.EQ, new NumericCompound<Double>(
            new FunctionExpression<>(MathFunctions.SIN, x),
            NumericOperator.DIV, new FunctionExpression<>(
                MathFunctions.COS, x))));
    TstUtil.runTest(expr, Result.SAT);
  }
  
  
}
