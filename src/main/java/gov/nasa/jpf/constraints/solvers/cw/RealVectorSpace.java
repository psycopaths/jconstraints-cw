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

import gov.nasa.jpf.constraints.api.Variable;
import gov.nasa.jpf.constraints.types.BuiltinTypes.DoubleType;
import gov.nasa.jpf.constraints.types.IntegerType;
import gov.nasa.jpf.constraints.types.RealType;
import gov.nasa.jpf.constraints.types.Type;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Real vector space that uses {@link SymbolicInteger}s and {@link SymbolicReal}s
 * to label its dimensions.
 * 
 * The RealVectorSpace maps the dimension labels to unique integer indices,
 * which allows {@link RealVector}s to store their values in arrays (instead of
 * maps).
 * 
 * @author Peter Dinges <pdinges@acm.org>
 */
class RealVectorSpace {

  private final List<Variable<?>> dimensions;
  private final Map<Object, Integer> dimensionIndices;

  private RealVectorSpace(List<Variable<?>> dimensions) {
    this.dimensions = Collections.unmodifiableList(dimensions);
    this.dimensionIndices = new HashMap<>(dimensions.size());
    for (int i = 0; i < dimensions.size(); ++i) {
      dimensionIndices.put(dimensions.get(i), i);
    }
  }

  public List<Variable<?>> dimensions() {
    return dimensions;
  }

  public int indexOf(Object dimension) {
    Integer index = dimensionIndices.get(dimension);
    if (index == null) {
      throw new IllegalArgumentException("Unknown dimension: " + dimension);
    }
    return index;
  }

  public RealVector.Builder makeVector() {
    return new RealVector.Builder(this);
  }

  public RealVector.Builder makeVector(RealVector initialValues) {
    return new RealVector.Builder(initialValues);
  }

  public MutableRealVector makeMutableVector() {
    return new MutableRealVector(this);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder(RealVectorSpace.class.getSimpleName());
    sb.append("{");
    for (int i = 0; i < dimensions.size(); ++i) {
      if (i > 0) {
        sb.append(", ");
      }
      sb.append(labelDimension(dimensions.get(i)));
    }
    sb.append("}");
    return sb.toString();
  }

  public static String labelDimension(Variable<?> dimension) {
    if (dimension.getType() instanceof IntegerType<?> ||
        dimension.getType() instanceof RealType<?>) {
      return dimension.getName();
    } else {
      throw new IllegalArgumentException("Unknown dimension type: " + dimension);
    }
  }

  public static RealVectorSpace forDimensions(Set<Variable<?>> variables) {
    for (Variable<?> var : variables) {
      Type<?> t = var.getType();
      if (!(t instanceof IntegerType<?> || t instanceof RealType<?>)) {
        throw new IllegalArgumentException("Variable " + var.getName() + " has unsupported variable type: " + t);
      }
    }
    List<Variable<?>> sortedDimensions = new ArrayList<>(variables);
    Collections.sort(sortedDimensions, ALPHABETICAL);
    return new RealVectorSpace(sortedDimensions);
  }

  public static RealVectorSpace forDimensions(Variable<?>... variables) {
    return forDimensions(new HashSet<Variable<?>>(Arrays.asList(variables)));
  }

  public static RealVectorSpace extend(RealVectorSpace space, Set<Variable<?>> extraVariables) {
    for (Object var : space.dimensions) {
      if (extraVariables.contains(var)) {
        throw new IllegalArgumentException("Extra variable already in space: " + var);
      }
    }
    List<Variable<?>> sortedDimensions = new ArrayList<>(space.dimensions);
    sortedDimensions.addAll(extraVariables);
    Collections.sort(sortedDimensions, ALPHABETICAL);
    return new RealVectorSpace(sortedDimensions);
  }

  private final static Comparator<Object> ALPHABETICAL = new Comparator<Object>() {
    @Override
    public int compare(Object arg0, Object arg1) {
      return arg0.toString().compareTo(arg1.toString());
    }
  };
}
