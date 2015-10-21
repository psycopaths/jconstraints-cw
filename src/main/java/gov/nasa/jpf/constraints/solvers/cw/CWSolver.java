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
//
// Copyright (C) 2014 The University of Illinois Board of Trustees.
// All Rights Reserved.
//
// The software in this Java package is distributed under the MIT License.
// See the MIT-LICENSE file at the top of the distribution directory tree
// for the complete license.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
package gov.nasa.jpf.constraints.solvers.cw;

import java.lang.annotation.RetentionPolicy;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Random;
import java.util.Set;
import java.util.logging.Logger;

import org.apache.commons.math3.random.MersenneTwister;

import gov.nasa.jpf.constraints.api.ConstraintSolver.Result;
import gov.nasa.jpf.constraints.api.ConstraintSolver;
import gov.nasa.jpf.constraints.api.Expression;
import gov.nasa.jpf.constraints.api.SolverContext;
import gov.nasa.jpf.constraints.api.Valuation;
import gov.nasa.jpf.constraints.api.ValuationEntry;
import gov.nasa.jpf.constraints.api.Variable;
import gov.nasa.jpf.constraints.expressions.Constant;
import gov.nasa.jpf.constraints.expressions.LogicalOperator;
import gov.nasa.jpf.constraints.expressions.NumericBooleanExpression;
import gov.nasa.jpf.constraints.expressions.NumericComparator;
import gov.nasa.jpf.constraints.expressions.NumericCompound;
import gov.nasa.jpf.constraints.expressions.NumericOperator;
import gov.nasa.jpf.constraints.expressions.PropositionalCompound;
import gov.nasa.jpf.constraints.solvers.ConstraintSolverFactory;
import gov.nasa.jpf.constraints.solvers.cw.CWConfig;
import gov.nasa.jpf.constraints.solvers.cw.CWSolverContext;
import gov.nasa.jpf.constraints.solvers.cw.CWSolver.CWSolverBuilder;
import gov.nasa.jpf.constraints.types.BuiltinTypes;
import gov.nasa.jpf.constraints.types.ConcreteIntegerType;
import gov.nasa.jpf.constraints.types.IntegerType;
import gov.nasa.jpf.constraints.types.RealType;
import gov.nasa.jpf.constraints.types.Type;
import gov.nasa.jpf.constraints.types.BuiltinTypes.BoolType;
import gov.nasa.jpf.constraints.types.BuiltinTypes.DoubleType;
import gov.nasa.jpf.constraints.types.BuiltinTypes.SInt16Type;
import gov.nasa.jpf.constraints.util.ExpressionClassifier;
import gov.nasa.jpf.constraints.util.ExpressionUtil;

public class CWSolver extends ConstraintSolver {
  public static class CWSolverBuilder {
    private Properties config;
    
    private boolean seedSet = false;
    private long seed;
    private String linearSolverStr = CWConfig.LINEAR_SOLVER.getDefaultValue();
    private int iterationsPerConstraint = CWConfig.ITERATIONS_PER_CONSTRAINT.getDefaultValue();
    private int neighborsGeneratedPerIteration = CWConfig.NEIGHBORS_GENERATED_PER_ITERATION.getDefaultValue();
    private float tabuIterationsPerVariable = CWConfig.TABU_ITERATIONS_PER_VARIABLE.getDefaultValue();
    private int minTabuIterations = CWConfig.MIN_TABU_ITERATIONS.getDefaultValue();
    private boolean enableSeeding = CWConfig.ENABLE_SEEDING.getDefaultValue();
    private boolean enableBisection = CWConfig.ENABLE_BISECTION.getDefaultValue();
    private int randomizationRadiusExp = CWConfig.RANDOMIZATION_RADIUS_EXPONENT.getDefaultValue();
    
    public CWSolverBuilder setLinearSolver(String linearSolver) {
      this.linearSolverStr = linearSolver;
      return this;
    }

    public CWSolverBuilder setIterationsPerConstraint(int iterationsPerConstraint) {
      this.iterationsPerConstraint = iterationsPerConstraint;
      return this;
    }

    public CWSolverBuilder setNeighborsGeneratedPerIteration(int neighborsGeneratedPerIteration) {
      this.neighborsGeneratedPerIteration = neighborsGeneratedPerIteration;
      return this;
    }

    public CWSolverBuilder setTabuIterationsPerVariable(float tabuIterationsPerVariable) {
      this.tabuIterationsPerVariable = tabuIterationsPerVariable;
      return this;
    }

    public CWSolverBuilder setMinTabuIterations(int minTabuIterations) {
      this.minTabuIterations = minTabuIterations;
      return this;
    }

    public CWSolverBuilder setEnableSeeding(boolean enableSeeding) {
      this.enableSeeding = enableSeeding;
      return this;
    }

    public CWSolverBuilder setEnableBisection(boolean enableBisection) {
      this.enableBisection = enableBisection;
      return this;
    }

    public CWSolverBuilder setRandomizationRadiusExp(int randomizationRadiusExp) {
      this.randomizationRadiusExp = randomizationRadiusExp;
      return this;
    }    
    
    public CWSolverBuilder setConfig(Properties conf) {
      this.config = conf;
      return this;
    }  
    
    public CWSolverBuilder setSeed(long seed) {
      this.seed = seed;
      this.seedSet = true;
      return this;
    }  
    
    public CWSolver buildSolver() {
      ConstraintSolverFactory csFac = (this.config != null) ? new ConstraintSolverFactory(config) : new ConstraintSolverFactory();
      ConstraintSolver linearSolver = csFac.createSolver(linearSolverStr);
      if(!this.seedSet)
        return new CWSolver(linearSolver,
                            this.iterationsPerConstraint,
                            this.neighborsGeneratedPerIteration,
                            this.tabuIterationsPerVariable,
                            this.minTabuIterations,
                            this.enableSeeding,
                            this.enableBisection,
                            this.randomizationRadiusExp);
      else
        return new CWSolver(linearSolver,
            this.iterationsPerConstraint,
            this.neighborsGeneratedPerIteration,
            this.tabuIterationsPerVariable,
            this.minTabuIterations,
            this.enableSeeding,
            this.enableBisection,
            this.randomizationRadiusExp,
            this.seed);
    }
  }
  
  
  /**
   * Produces a lot of log data atm
   */
  private final boolean USE_LOGGER = false;
  
  protected ConstraintSolver linearSolver;
  /**
   * Number of algorithm iterations added for each constraint in the path
   * condition. The paper speaks of <em>steps</em> instead of iterations and
   * denotes this constant as <em>I</em>. 
   */
  private final int ITERATIONS_PER_CONSTRAINT;
  /**
   * Number of random neighbor {@link RealVector points} to generate per iteration.
   * The paper calls points <em>environments</em> and denotes this constant
   * as <em>R</em>.  
   */
  private final int NEIGHBORS_GENERATED_PER_ITERATION;
  /**
   * Multiplier for computing how many iterations a variable is tabu if changing
   * it failed to yield a better point. The paper uses this constant and
   * {@link #MIN_TABU_ITERATIONS} to compute the number of tabu iterations
   * <em>T</em>.
   */
  private final float TABU_ITERATIONS_PER_VARIABLE;
  /**
   * Minimum number of iterations that a variable is tabu if changing it failed
   * to yield a better point. See {@link #TABU_ITERATIONS_PER_VARIABLE} above.
   */
  private final int MIN_TABU_ITERATIONS;
  /**
   * Enable trying the hard-coded {@link #SEEDED_VALUES} for the changed variable?
   */
  private final boolean ENABLE_SEEDING;
  /**
   * Enable estimating additional neighbors that <em>should</em> satisfy one
   * of the constraints violated by the current variable if the constraint was
   * linear?  See the discussion of Finding Neighbors in section 3.3 of the paper.
   */
  public final boolean ENABLE_BISECTION;
  
  private final int RANDOMIZATION_RADIUS_EXPONENT;
  /** Error value used if the error function of a constraint returns NaN */
  private final static double NAN_ERROR = 0.125 * Double.MAX_VALUE;
  private final static double[] SEEDED_VALUES =
    { 0.0, 1.0, -1.0, Double.MAX_VALUE, Double.MIN_VALUE };

  private final Random random;
  
  private Logger logger = Logger.getLogger("constraints");
  
  public CWSolver(ConstraintSolver linearSolver, 
      int iterationsPerConstraint, 
      int neighborsGeneratedPerIteration, 
      float tabuIterationsPerVariable,
      int minTabuIterations,
      boolean enableSeeding, 
      boolean enableBisection, 
      int randomizationRadiusExp,
      long seed) {
    this.linearSolver = linearSolver;
    this.ITERATIONS_PER_CONSTRAINT = iterationsPerConstraint;
    this.NEIGHBORS_GENERATED_PER_ITERATION = neighborsGeneratedPerIteration;
    this.TABU_ITERATIONS_PER_VARIABLE = tabuIterationsPerVariable;
    this.MIN_TABU_ITERATIONS = minTabuIterations;
    this.ENABLE_BISECTION = enableBisection;
    this.ENABLE_SEEDING = enableSeeding;
    this.RANDOMIZATION_RADIUS_EXPONENT = randomizationRadiusExp;
    this.random = new Random(seed);
  }
  public CWSolver(ConstraintSolver linearSolver, 
                  int iterationsPerConstraint, 
                  int neighborsGeneratedPerIteration, 
                  float tabuIterationsPerVariable,
                  int minTabuIterations,
                  boolean enableSeeding, 
                  boolean enableBisection, 
                  int randomizationRadiusExp) {
    this(linearSolver,iterationsPerConstraint,
        neighborsGeneratedPerIteration, tabuIterationsPerVariable, 
        minTabuIterations, enableSeeding, enableBisection, randomizationRadiusExp, new MersenneTwister().nextLong());
  }
  
  @Override
  public SolverContext createContext() {
    return new CWSolverContext(this);
  }
  
  
  @Override
  public Result solve(Expression<Boolean> constraint, Valuation result) {
    
    List<Expression<Boolean>> linear = new ArrayList<>();
    List<Expression<Boolean>> nonlinear = new ArrayList<>();
    
    ExpressionClassifier.splitConstraintsInto(constraint, linear, nonlinear);
    
    Expression<Boolean> linearConstraint = (linear.size() > 0) ? ExpressionUtil.and(linear) : null;
    if(linearConstraint == null) { //TODO: anything else we can do here?
      return Result.DONT_KNOW;
    }
    
    
    //if there are no linear constraints, the best we can do is just to use an arbitrary valuation?
    //Expression<Boolean> linearConstraint =  ExpressionUtil.and(linear);

    Expression<Boolean> nonlinearConstraint = (nonlinear.size() > 0) ? ExpressionUtil.and(nonlinear) : null;
    //Expression<Boolean> nonlinearConstraint = ExpressionUtil.and(nonlinear);
    Valuation linearValuation = new Valuation();
    Result linSol = this.linearSolver.solve(linearConstraint, linearValuation);
    
    if(linSol != Result.SAT) { // lines 4--5
      return linSol;
    }

    if(nonlinearConstraint == null) { //no nonlinear constraints are present
      for(ValuationEntry<?> entry : linearValuation)
        result.addEntry(entry);
      return Result.SAT;
    }
    
    return solveWithAdaptiveVariableSearch(linearConstraint, linearValuation, nonlinearConstraint, result);
  }
  
  private static <T> Set<T> union(Collection<T> collection1, Collection<T> collection2) {
    Set<T> result = new HashSet<>();
    result.addAll(collection1);
    result.addAll(collection2);
    return result;
  }
  
  
  //TODO: note that adaptive variable search can ONLY return SAT or DONT_KNOW atm
  private Result solveWithAdaptiveVariableSearch(Expression<Boolean> linearConstraint, Valuation linearSolution, Expression<Boolean> nonlinearConstraint, Valuation result) {

    // Partially based on the paper "Yet Another Local Search Method for Constraint Solving"
    // by Philippe Codognet and Daniel Diaz, 2001.

    // NOTE: End-of-line comments give the variable names and line numbers
    //       used in Algorithms 1 and 2 of the paper.
    
    // Interpret every variable in the path condition as a dimension
    // in a real-valued vector space.
    Set<Variable<?>> vars = union(ExpressionUtil.freeVariables(linearConstraint), ExpressionUtil.freeVariables(nonlinearConstraint));
    RealVectorSpace vectorSpace =
        RealVectorSpace.forDimensions(vars);
    int numberOfVariables = vectorSpace.dimensions().size();

    RealVector p = makeVectorFromSolutions(vectorSpace, linearSolution);  // Alg. 1: \alpha, line 4

    printDebug(CWSolver.class, "Linear PC: ", linearConstraint);
    printDebug(CWSolver.class, "Linear PC solution: ", p);
    printDebug(CWSolver.class, "Solving non-linear PC\n", nonlinearConstraint);

    
    List<Expression<Boolean>> nonLinearConstraintsLst = ExpressionClassifier.splitToConjuncts(nonlinearConstraint);
    
    // Initialize lookup tables and local variables
    Map<Expression<Boolean>, BitSet> variableIndicesByConstraint = new IdentityHashMap<>(nonLinearConstraintsLst.size());
    @SuppressWarnings("unchecked")
    List<Expression<Boolean>>[] constraintsByVariableIndex = new List[numberOfVariables];
    populateLookupTables(vectorSpace, nonLinearConstraintsLst, constraintsByVariableIndex, variableIndicesByConstraint);

    double[] errorByVariables = new double[numberOfVariables];  // Alg. 1: \epsilon
    int[] tabuVariables = new int[numberOfVariables];  // Alg. 1: \tau

    int iterationCount = 1;  // i
    int iterationLimit = nonLinearConstraintsLst.size() * ITERATIONS_PER_CONSTRAINT;  // Alg. 1: I

    
    //iterate as long as non linear constraint is not satisfied
    while (!nonlinearConstraint.evaluate(p.convertToJconstraintsValuation())) {  // Alg. 1: line 7
      if (iterationCount > iterationLimit) {  // Alg. 1: line 8
        printDebug(CWSolver.class, "Could not find solution within ", iterationLimit, " iterations");
        return Result.DONT_KNOW;
      }
      ++iterationCount;  // Alg. 1: line 9

      // Compute errors
      double errorAtP = 0.0;  // Alg. 1: e_\alpha
      Arrays.fill(errorByVariables, 0.0);  // Alg. 1: line 14
      for (Expression<Boolean> c : nonLinearConstraintsLst) {  // Alg. 1: lines 15--20
        if(c instanceof PropositionalCompound) {
          System.out.println("propositional compound. Skipping");
          continue;
        }
          
        //TODO: fix the list; it must be composed of NumericBooleanExpressions (stronger type)
        if(!(c instanceof NumericBooleanExpression))
          throw new IllegalStateException("constraint must be " + NumericBooleanExpression.class.getName() + " and not of type " + c.getClass().getName());
        
        NumericBooleanExpression nc = (NumericBooleanExpression)c;
        
        double e = computeError(nc, p.convertToJconstraintsValuation());
        errorAtP += e;
        incrementElements(errorByVariables, variableIndicesByConstraint.get(c), e);
      }
      printDebug(CWSolver.class, "p = ", p, " -> error ", errorAtP);

      // Try to find a better solution by modifying the "worst" non-tabu variable
      int wiggleVarIndex = indexOfMaxIgnoringTabu(errorByVariables, tabuVariables);
      if (wiggleVarIndex == -1) { // All variables might be tabu,  Alg. 1: lines 10--13
        for (int i = 0; i < tabuVariables.length; ++i) {
          p = makeRandomNeighborInPolytope(p, linearConstraint, vectorSpace.dimensions().get(i));
          if (p == null) { //no point could be found, so dont_know?
            return Result.DONT_KNOW;
          }
          tabuVariables[i] = 0;
        }
        printDebug(CWSolver.class, "All variables are tabu.  Took random step to ", p);
        continue;
      }
      Variable<?> wiggleVar = vectorSpace.dimensions().get(wiggleVarIndex);  // Alg. 1: x, line 21
      printDebug(CWSolver.class, "Wiggling ", wiggleVar);

      // Find best neighbor (Algorithm 3)
      
      double minError = Double.POSITIVE_INFINITY;  // Alg. 3: e_\mu
      RealVector minNeighbor = null;
      for (int i = 0; i < NEIGHBORS_GENERATED_PER_ITERATION; ++i) {
        RealVector q = makeRandomNeighborInPolytope(p, linearConstraint, wiggleVar);  // Alg. 3: \beta
        RealVector r = null;  // Alg. 3: \gamma

        if (q == null) { // No random neighbor could be found
          break;
        }

        double errorAtQ = computeError(nonLinearConstraintsLst, q.convertToJconstraintsValuation());  // Alg. 3: e_\beta, line 7
        double errorAtR = Double.POSITIVE_INFINITY;  // Alg. 3: e_\gamma

        if (ENABLE_BISECTION && constraintsByVariableIndex[wiggleVarIndex] != null) {  // Alg. 3: line 5
          // Pick a random unsatisfied constraint
          List<Expression<Boolean>> constraintsForVar = new ArrayList<>(constraintsByVariableIndex[wiggleVarIndex]);
          Collections.shuffle(constraintsForVar);
          Expression<Boolean> constraint = null;
          for (int k = 0; k < constraintsForVar.size(); ++k) {
            constraint = constraintsForVar.get(k);
            
            boolean sat = constraint.evaluate(p.convertToJconstraintsValuation());
            if (!sat) {
              break; 
            }
          }
          Number valueAtP = evaluateAndSubtractSides((NumericBooleanExpression)constraint, p.convertToJconstraintsValuation());
          Number valueAtQ = evaluateAndSubtractSides((NumericBooleanExpression)constraint, q.convertToJconstraintsValuation());
          r = linearlyEstimateZero(p, valueAtP, q, valueAtQ, ((NumericBooleanExpression)constraint).getComparator());  // Alg. 3: line 6

          boolean sat = linearConstraint.evaluate(r.convertToJconstraintsValuation());
          
          if (sat) {
            errorAtR = computeError(nonLinearConstraintsLst, r.convertToJconstraintsValuation());
          }
        }

        printDebug(CWSolver.class, "Random neighbors");
        printDebug(CWSolver.class, "    q = ", q, " -> error ", errorAtQ);
        printDebug(CWSolver.class, "    r = ", r, " -> error ", errorAtR);

        if (errorAtQ < minError) {  // Alg. 3: lines 9--12
          minError = errorAtQ;
          minNeighbor = q;
        }
        if (errorAtR < minError) { // Alg. 3: lines 13--16
          minError = errorAtR;
          minNeighbor = r;
        }
      }

      if (ENABLE_SEEDING) {
        for (double seed : SEEDED_VALUES) {
          RealVector s = vectorSpace.makeVector(p).set(wiggleVar, seed).build();
          double errorAtS = computeError(nonLinearConstraintsLst, s.convertToJconstraintsValuation());
          
          boolean sat = linearConstraint.evaluate(s.convertToJconstraintsValuation());
          if (sat && errorAtS < minError) {
            minError = errorAtS;
            minNeighbor = s;
          }
        }
      }

      if (minError < errorAtP) {  // Alg. 1: lines 23--27
        printDebug(CWSolver.class, "Found new neighbor");
        p = minNeighbor;
        decrementElements(tabuVariables);
      } else {  // Alg 1: lines 27--29
        printDebug(CWSolver.class, "Could not find better neighbor");
        tabuVariables[wiggleVarIndex] = Math.max(Math.round(TABU_ITERATIONS_PER_VARIABLE * vectorSpace.dimensions().size()), MIN_TABU_ITERATIONS);
        printDebug(CWSolver.class, "Tabu ", Arrays.toString(tabuVariables));
      }
    }
    printDebug(CWSolver.class, "Found solution: ", p);
    for(Variable<?> v : vars) {
      double vectorVal = p.get(v);
      if(v.getType() instanceof IntegerType<?>) {
        result.setValue((Variable<Integer>)v, (int)vectorVal);
      } else if(v.getType() instanceof RealType<?>) {
        result.setValue((Variable<Double>)v, vectorVal);
      }
    }
    System.out.println("cwsolver solution: " + result.toString());
    return Result.SAT;
  }
  
  private static void typeCheck(Variable<?> var) {
    Type<?> t = var.getType();
    if(!(t instanceof IntegerType<?> || t instanceof RealType<?>))
      throw new AssertionError("Unknown variable type: " + var);
  }
  
  /**
   * Returns a point in the {@code vectorSpace} that is described by the
   * solution values attached to the given set of {@code variables}.
   */
  private static RealVector makeVectorFromSolutions(RealVectorSpace vectorSpace, Valuation valuation) {
    RealVector.Builder builder = vectorSpace.makeVector();
    
    if (valuation.entries().size() == 0) {
      return builder.build();
    }
    
    for(ValuationEntry<?> sol : valuation) {
      Variable<?> var = sol.getVariable();
      typeCheck(var);
      Object value = sol.getValue();
      if(value == null)
        throw new AssertionError("Variable " + var + " has no attached solution");
      
      if(!(value instanceof Number))
        throw new AssertionError("Unsupported value type" + value.getClass().getName());
      
      builder.set(var, ((Number)value).doubleValue());
    }
    return builder.build();
  }
  
  static void populateLookupTables(RealVectorSpace space, List<Expression<Boolean>> constraints, List<Expression<Boolean>>[] constraintsByVariableIndex, Map<Expression<Boolean>, BitSet> variableIndicesByConstraint) {
    int numberOfVariables = space.dimensions().size();
    for (Expression<Boolean> constraint : constraints) {
      BitSet variableIndices = new BitSet(numberOfVariables);
      variableIndicesByConstraint.put(constraint, variableIndices);
      for (Variable<?> variable : ExpressionUtil.freeVariables(constraint)) {
        int varIndex = space.indexOf(variable);
        variableIndices.set(varIndex);
        if (constraintsByVariableIndex[varIndex] == null) {
          constraintsByVariableIndex[varIndex] = new ArrayList<>();
        }
        constraintsByVariableIndex[varIndex].add(constraint);
      }
    }
  }
  
  /**
   * Returns the cumulative error score for the constraints in the given path
   * condition at the given valuation point.
   * 
   * <p>
   * This is Algorithm 2 in the paper.
   */
  static double computeError(List<Expression<Boolean>> constraints, Valuation valuation) {
    double totalError = 0;
    for (Expression<Boolean> c : constraints) {
      if(c instanceof PropositionalCompound) {
        System.out.println("propositional compound. Skipping");
        continue;
      }
      if(!(c instanceof NumericBooleanExpression)) {
        throw new IllegalStateException("must be of type " + NumericBooleanExpression.class.getName() + "; not " + c.getClass().getName());
      }
      totalError += computeError((NumericBooleanExpression)c, valuation);
    }
    return totalError;
  }

  
  //DELETE?
  static double computeError(Expression<Boolean> constraint, Valuation valuation) {
    NumericBooleanExpression bExp = new NumericBooleanExpression(constraint, NumericComparator.EQ, new Constant<>(BuiltinTypes.BOOL, true));
    return computeError(bExp.getComparator(), evaluateAndSubtractSides(bExp, valuation));
  }
  
  static double computeError(NumericBooleanExpression constraint, Valuation valuation) {
    NumericBooleanExpression bExp = (NumericBooleanExpression)constraint;
    return computeError(bExp.getComparator(), evaluateAndSubtractSides(bExp, valuation));
  }

  static double computeError(NumericComparator comparator, Number value) {
    return computeError(comparator, value.doubleValue());
  }

  static double computeError(NumericComparator comparator, double value) {
    if (Double.isNaN(value)) {
      return NAN_ERROR;
    }
    switch (comparator) {
    case EQ:
      return Math.abs(value);
    case GE:
      return value >= 0 ? 0 : -value + 1;
    case GT:
      return value > 0 ? 0 : -value + 1;
    case LE:
      return value <= 0 ? 0 : value + 1;
    case LT:
      return value < 0 ? 0 : value + 1;
    case NE:
      return value != 0 ? 0 : 1;
    }
    throw new AssertionError();
  }
  
  public static Number evaluateAndSubtractSides(NumericBooleanExpression constraint, Valuation valuation) {
    Object leftVal = constraint.getLeft().evaluate(valuation);
    Object rightVal = constraint.getRight().evaluate(valuation);

    if(!(leftVal instanceof Number && rightVal instanceof Number) &&
       !(leftVal instanceof Boolean && rightVal instanceof Boolean))
      throw new IllegalStateException("left and right hand side must evaluate to a number or boolean; not Ltype: " + leftVal.getClass().getName() + " RType: " + rightVal.getClass().getName());
    Number l; //DELETE?
    if(leftVal instanceof Boolean) {
      Boolean lv = (Boolean)leftVal;
      l = (lv) ? 1.0 : 0.0;
    } else
      l = (Number)leftVal;
    boolean leftIsNegative = l.doubleValue() < 0;
    
    Number r; //DELETE?
    if(rightVal instanceof Boolean) {
      Boolean rv = (Boolean)rightVal;
      r = (rv) ? 1.0 : 0.0;
    } else
      r = (Number)rightVal;
    boolean rightIsNegative = r.doubleValue() < 0;

    //TODO: not sure if this is correct.. double check with original implementation
    Number result = l.doubleValue() - r.doubleValue();
    if (leftIsNegative == rightIsNegative) {
      // Can neither overflow nor underflow
      return result;
    }

    boolean resultIsNegative = result.doubleValue() < 0;
    boolean underflowOccurred = leftIsNegative && !resultIsNegative;
    boolean overflowOccurred = rightIsNegative && resultIsNegative;
    if (!underflowOccurred && !overflowOccurred) {
      return result;
    }

    if (result instanceof Integer) {
      return (underflowOccurred ? Integer.MIN_VALUE : Integer.MAX_VALUE);
    } else if (result instanceof Double) {
      return (underflowOccurred ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY);
    } else {
      throw new AssertionError("Result has unknown number type: " + result);
    }
  }
  
  static void incrementElements(double[] values, BitSet selectedIndices, double increment) {
    for (int i = selectedIndices.nextSetBit(0); i >= 0; i = selectedIndices.nextSetBit(i+1)) {
      values[i] += increment;
    }
  }

  static void decrementElements(int[] values) {
    for (int i = 0; i < values.length; ++i) {
      if (values[i] > 0) {
        --values[i];
      }
    }
  }
  static int indexOfMaxIgnoringTabu(double[] values, int[] tabuIndices) {
    double currentMax = Double.NEGATIVE_INFINITY;
    int currentMaxIndex = -1;

    for (int i = 0; i < values.length; ++i) {
      if (tabuIndices[i] > 0) {
        continue;
      }
      if (currentMax < values[i]) {
        currentMaxIndex = i;
        currentMax = values[i];
      }
    }
    return currentMaxIndex;
  }
  
  /**
   * Returns a new point in the {@code polytope} that differs from the given
   * {@code point} in the given {@code variable} (= dimension), or {@code null}
   * if no such point could be found.
   * 
   * @see Algorithm 4 in the paper.
   */
  private RealVector makeRandomNeighborInPolytope(RealVector point, Expression<Boolean> polytope, Variable<?> variable) {
    // Implement this as Dikin walk, see Kannan and Narayanan "Random Walks
    // on Polytopes and an Affine Interior Point Method for Linear Programing".
    RealVector.Builder deltaBuilder = point.space().makeVector();
    RealVector neighbor;

    final int e = 2 * RANDOMIZATION_RADIUS_EXPONENT;
    for (int i = 0; i < 25; ++i) {
      double g = random.nextDouble();
      int j = random.nextInt(e) - RANDOMIZATION_RADIUS_EXPONENT;
      double delta = (j > 0) ? g * (1 << j) : 1.0 / (g * (1 << -j));
      deltaBuilder.set(variable, delta);
      neighbor = RealVector.plus(point, deltaBuilder.build());

      boolean satisfied = polytope.evaluate(neighbor.convertToJconstraintsValuation());
      
      if (satisfied) {
        return neighbor;
      }
    }
    printDebug(CWSolver.class, "Failed to find random neighbor in polytope");
    return null;
  }
  
  /**
   * Returns a point {@code r} where the linear function {@code f} with
   * {@code f(p) = valueAtP} and {@code f(q) = valueAtQ} is zero, that is,
   * {@code f(r) = 0}.
   * 
   * <p> 
   * NOTE: The returned point does not necessarily lie within the polytope.
   * 
   * <p>
   * This is the BisectionStep function that appears in the paper.
   */
  private RealVector linearlyEstimateZero(RealVector p, Number valueAtP, RealVector q, Number valueAtQ, NumericComparator comparator) {
    double vP = valueAtP.doubleValue();
    double vQ = valueAtQ.doubleValue();
    double t0 = (vP != vQ) ? -vP / (vQ - vP) : random.nextGaussian();
    if (comparator != NumericComparator.EQ) {
      t0 *= (1.0 + (random.nextInt(250) * 0.01));
      if (t0 == 0) {
        t0 = random.nextGaussian();
      }
    }
    return RealVector.plus(RealVector.times(t0, RealVector.minus(q, p)), p);
  }
  
  /*************************************************************************
   * 
   * 
   * Debugging
   * TODO: change to Logger.....
   * 
   */
  public void printDebug(Class<?> source, Object message) {
    printDebug(source, message.toString());
  }

  public void printDebug(Class<?> source, Object message1, Object message2) {
    printDebug(source, message1.toString() + message2);
  }

  public void printDebug(Class<?> source, Object message1, Object message2, Object message3) {
    printDebug(source, message1.toString() + message2 + message3);
  }

  public void printDebug(Class<?> source, Object message1, Object message2, Object message3, Object message4) {
    printDebug(source, message1.toString() + message2 + message3 + message4);
  }


  public void printDebug(Class<?> source, Object message1, Object message2, Object message3, Object message4, Object... messages) {
    StringBuffer sb = new StringBuffer();
    sb.append(message1).append(message2).append(message3).append(message4);
    for (Object m : messages) {
      sb.append(m);
    }
    printDebug(source, sb.toString());
  }

  public void printDebug(Object source, String message) {
    printDebug(source.getClass(), message);
  }

  public void printDebug(Class<?> klass, String message) {
    if(USE_LOGGER)
      logger.info("[CW] - " + klass.getSimpleName() + ": " + message);
  }

  public static String defaultStringRepresentation(Object this_, Object... fieldValues) {
    StringBuilder sb = new StringBuilder(this_.getClass().getSimpleName());
    sb.append("{");
    for (int i = 0; i < fieldValues.length; ++i) {
      if (i > 0) {
        sb.append(", ");
      }
      sb.append(fieldValues[i]);
    }
    sb.append("}");
    return sb.toString();
  }
}
