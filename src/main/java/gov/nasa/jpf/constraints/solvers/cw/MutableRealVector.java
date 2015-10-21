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

import gov.nasa.jpf.constraints.api.Variable;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

class MutableRealVector {

  private final RealVectorSpace space;
  private final double[] values;

  public MutableRealVector(RealVectorSpace space) {
    this.space = space;
    this.values = new double[space.dimensions().size()];
  }

  public MutableRealVector(RealVectorSpace space, double defaultValue) {
    this(space);
    Arrays.fill(this.values, defaultValue);
  }

  public MutableRealVector(RealVectorSpace space, double[] values) {
    assert space.dimensions().size() == values.length;
    this.space = space;
    this.values = Arrays.copyOf(values, values.length);  // Defensive copy
  }

  public RealVectorSpace space() {
    return space;
  }

  public double get(Object dimension) {
    return values[space.indexOf(dimension)];
  }

  public MutableRealVector set(Object dimension, double value) {
    values[space.indexOf(dimension)] = value;
    return this;
  }

  public boolean isZeroExcept(Set<Object> dimensions) {
    for (int i = 0; i < values.length; ++i) {
      if (values[i] != 0
          && !(dimensions.contains(space.dimensions().get(i)))) {
        return false;
      }
    }
    return true;
  }

  // Operations

  public MutableRealVector clear() {
    for (int i = 0; i < values.length; ++i) {
      this.values[i] = 0.0;
    }
    return this;
  }

  public MutableRealVector add(MutableRealVector other) {
    checkSameSpace(this, other);
    int n = this.values.length;
    for (int i = 0; i < n; ++i) {
      this.values[i] += other.values[i];
    }
    return this;
  }

  public MutableRealVector subtract(MutableRealVector other) {
    checkSameSpace(this, other);
    int n = this.values.length;
    for (int i = 0; i < n; ++i) {
      this.values[i] -= other.values[i];
    }
    return this;
  }

  public MutableRealVector multiply(double scalar) {
    for (int i = 0; i < values.length; ++i) {
      values[i] *= scalar;
    }
    return this;
  }

  public MutableRealVector divide(double scalar) {
    for (int i = 0; i < values.length; ++i) {
      values[i] /= scalar;
    }
    return this;
  }

  @Override
  public boolean equals(Object other) {
    if (this == other) {
      return true;
    }
    if (!(other instanceof MutableRealVector)) {
      return false;
    }
    MutableRealVector otherVector = (MutableRealVector) other;
    return this.space == otherVector.space
        && Arrays.equals(this.values, otherVector.values);
  }

  @Override
  public int hashCode() {
    int hash = space.hashCode() * 31;
    for (double val : values) {
      hash = hash * 31 + ((int) val);
    }
    return hash;
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder(MutableRealVector.class.getSimpleName());
    sb.append("{");
    List<Variable<?>> dimensions = space.dimensions();
    for (int i = 0; i < dimensions.size(); ++i) {
      if (i > 0) {
        sb.append(", ");
      }
      sb.append(RealVectorSpace.labelDimension(dimensions.get(i))).append(": ");
      sb.append(values[i]);
    }
    sb.append("}");
    return sb.toString();
  }

  private static void checkSameSpace(MutableRealVector a, MutableRealVector b) {
    if (a.space != b.space) {
      throw new IllegalArgumentException("Vectors must belong to the same space");
    }
  }
}
