/*******************************************************/
/* Copyright (c) 2015 by Artelys                       */
/* All Rights Reserved                                 */
/*******************************************************/
#pragma once

#include <vector>
#include "knitro.h"

namespace knitro {

class KTRISolver;

/**
 * Abstract base class for a new point callback. 
 * The CallbackFunction is called after each new estimate of the solution point.
 */
class KTRNewptCallback {
 public:

  /**
   * Virtual destructor required for virtual classes. 
   */
  virtual ~KTRNewptCallback() {
  }

  /**
   *
   * @param x
   * @param lambda
   * @param obj
   * @param c
   * @param objGrad
   * @param jac
   * @param solver
   * @return
   */
  virtual int CallbackFunction(const double* const x, const double* const lambda, double obj,
                               const double* const c, const double* const objGrad,
                               const double* const jac, KTRISolver * solver) = 0;

  /**
   * @param x
   * @param obj
   * @param objGrad
   * @param jac
   * @param solver
   * @return
   */ 
  virtual int CallbackFunctionLSQ(const double* const x, double obj, const double* const objGrad,
                                  const double* const jac, KTRISolver * solver) { return -1; }

};

}

