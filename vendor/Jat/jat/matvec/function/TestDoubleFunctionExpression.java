/* JAT: Java Astrodynamics Toolkit
 *
 * Copyright (c) 2003 United States Government as represented by the
 * administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 *
 * This file is part of JAT. JAT is free software; you can
 * redistribute it and/or modify it under the terms of the
 * NASA Open Source Agreement, version 1.3 or later.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * NASA Open Source Agreement for more details.
 *
 * You should have received a copy of the NASA Open Source Agreement
 * along with this program; if not, write to the NASA Goddard
 * Space Flight Center at opensource@gsfc.nasa.gov.
 */

package jat.matvec.function;

import jat.matvec.function.expressionParser.Evaluator;
import jat.matvec.data.Matrix;

public class TestDoubleFunctionExpression extends TestDoubleFunction {

  //private int argNumber;

  private String expression;
  private String[] variables;
  private Evaluator E = new Evaluator();

  public TestDoubleFunctionExpression(String exp,String[] vars) {
    argNumber = vars.length;
    setFunction(exp,vars);
  }

  public TestDoubleFunctionExpression(String exp,String vars) {
    argNumber = 1;
    String[] variable = new String[1];
    variable[0] = vars;
    setFunction(exp,variable);
  }

  private void setFunction(String exp,String[] vars) {
    expression = exp;
    variables = vars;
  }


  private void setVariableValues(double[] values) {
    for (int i = 0;i<variables.length;i++) {
      E.addVariable(variables[i],values[i]);
    }
    E.setExpression(expression);
  }

  private void setVariableValues(double value) {
    E.addVariable(variables[0],value);
    E.setExpression(expression);
  }

  private void setVariableValues(Matrix values) {
    for (int i = 0;i<variables.length;i++) {
      E.addVariable(variables[i],values.get(i,0));
    }
    E.setExpression(expression);
  }

  public boolean eval(double[] values) {
    checkArgNumber(values.length);
    setVariableValues(values);
    boolean bool = (((Double)(E.getValue())).doubleValue()>0.5) ? (true) : (false);
    return bool;
  }

  public boolean eval(double value) {
    checkArgNumber(1);
    setVariableValues(value);
    boolean bool = (((Double)(E.getValue())).doubleValue()>0.5) ? (true) : (false);
    return bool;
  }

  public boolean eval(Matrix values) {
    checkArgNumber(values.getRowDimension());
    setVariableValues(values);
    boolean bool = (((Double)(E.getValue())).doubleValue()>0.5) ? (true) : (false);
    return bool;
  }
}
