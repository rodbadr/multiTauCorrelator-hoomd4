#include "test_correlator_likh.h"
#include <math.h>
#include <iostream>

using namespace std;

/////////////////////////////////////////
// Correlator class
/////////////////////////////////////////
Correlator_Likh::Correlator_Likh(const unsigned int numcorrin,
                                 const unsigned int p_in,
                                 const unsigned int m_in) {
  setsize(numcorrin, p_in, m_in);
}

Correlator_Likh::~Correlator_Likh() {
  if (numcorrelators == 0) return;

  delete[] shiftA;
  delete[] shiftB;
  delete[] correlation;
  delete[] ncorrelation;
  delete[] accumulatorA;
  delete[] accumulatorB;
  delete[] naccumulator;
  delete[] insertindex;

  delete[] t;
  delete[] f;
}

void Correlator_Likh::setsize(const unsigned int numcorrin,
                              const unsigned int p_in, const unsigned int m_in) {

/*
  -numcorrin: The number of correlators (levels) in the multi-tau correlator hierarchy.
  -p_in: The number of time bins for each correlator.
  -m_in: The decimation factor (or averaging factor) between correlator levels. 
  After every m_in data points at one level, their average is passed to the next level.
*/
  numcorrelators = numcorrin;
  p = p_in;
  m = m_in;
  d_min = p / m;

  length = numcorrelators * p;

  shiftA = new double*[numcorrelators];
  shiftB = new double*[numcorrelators];
  correlation = new double*[numcorrelators];
  ncorrelation = new unsigned long int*[numcorrelators];
  accumulatorA = new double[numcorrelators];
  accumulatorB = new double[numcorrelators];
  naccumulator = new unsigned int[numcorrelators];
  insertindex = new unsigned int[numcorrelators];

  for (unsigned int j = 0; j < numcorrelators; ++j) {
    shiftA[j] = new double[p];
    shiftB[j] = new double[p];

    /* It can be optimized: Apart from correlator 0, correlation and
     * ncorrelation arrays only use p/2 values */
    correlation[j] = new double[p];
    ncorrelation[j] = new unsigned long int[p];
  }

  t = new double[length];
  f = new double[length];
}

void Correlator_Likh::initialize() {
  for (unsigned int j = 0; j < numcorrelators; ++j) {
    for (unsigned int i = 0; i < p; ++i) {
      shiftA[j][i] = -2E10;
      shiftB[j][i] = -2E10;
      correlation[j][i] = 0;
      ncorrelation[j][i] = 0;
    }
    accumulatorA[j] = 0.0;
    accumulatorB[j] = 0.0;
    naccumulator[j] = 0;
    insertindex[j] = 0;
  }

  for (unsigned int i = 0; i < length; ++i) {
    t[i] = 0;
    f[i] = 0;
  }

  npcorr = 0;
  kmax = 0;
  accvalA = 0;
  accvalB = 0;
}

void Correlator_Likh::add(const double wA, const double wB, const unsigned int k) {
/*
The add function is called every time you have a new value (e.g., at each timestep). It:

  -Inserts the new value into the correlator’s internal buffers.
  -Updates accumulators and, if enough values have been collected, recursively adds averaged values to higher-level correlators.
  -Updates the running sums needed for correlation calculation.
*/

  /// If we exceed the correlator side, the value is discarded
  if (k == numcorrelators) return;
  if (k > kmax) kmax = k;

  /// (1) Insert new value in shift array
  shiftA[k][insertindex[k]] = wA;
  shiftB[k][insertindex[k]] = wB;

  /// Add to average value
  if (k == 0) {
    accvalA += wA;
    accvalB += wB;
  }

  /// (3 & 4 Add to accumulator and, if needed, add to next correlator
  accumulatorA[k] += wA;
  accumulatorB[k] += wB;
  ++naccumulator[k];
  if (naccumulator[k] == m) { 
    double avgA = accumulatorA[k] / m;
    double avgB = accumulatorB[k] / m;
    add(avgA, avgB, k + 1);
    accumulatorA[k] = 0.0;
    accumulatorB[k] = 0.0;
    naccumulator[k] = 0;
  }

  /// Calculate correlation function
  unsigned int ind1 = insertindex[k];
  if (k == 0) {  /// First correlator is different
    int ind2 = ind1;
    for (unsigned int j = 0; j < p; ++j) {
      if (shiftA[k][ind2] > -1e10) {
        correlation[k][j] += shiftA[k][ind1] * shiftB[k][ind2];
        ++ncorrelation[k][j];
      }
      --ind2;
      if (ind2 < 0) ind2 += p;
    }
  } else {
    int ind2 = ind1 - d_min;
    for (unsigned int j = d_min; j < p; ++j) {
      if (ind2 < 0) ind2 += p;
      if (shiftA[k][ind2] > -1e10) {
        correlation[k][j] += shiftA[k][ind1] * shiftB[k][ind2];
        ++ncorrelation[k][j];
      }
      --ind2;
    }
  }

  ++insertindex[k];
  if (insertindex[k] == p) insertindex[k] = 0;
}

void Correlator_Likh::evaluate(const bool norm) {
  /*
  The evaluate function is called less frequently (e.g., at analysis intervals). It:

  -Processes all the accumulated data to compute the actual correlation function.
  -Fills the output arrays (t and f) with the time lags and corresponding correlation values.
  -Sets npcorr to the number of valid correlation points.

  \param norm: If norm is true, it normalizes the correlation values by the average value of the first correlator.
*/
  unsigned int im = 0;

  double aux = 0;
  if (norm && ncorrelation[0][0] > 0) {
    double meanA = accvalA / ncorrelation[0][0];
    double meanB = accvalB / ncorrelation[0][0];
    aux = meanA * meanB;
  }

  // First correlator
  for (unsigned int i = 0; i < p; ++i) {
    if (ncorrelation[0][i] > 0) {
      t[im] = i;
      f[im] = correlation[0][i] / ncorrelation[0][i] - aux;
      ++im;
    }
  }

  // Subsequent correlators
  for (int k = 1; k < kmax; ++k) {
    for (int i = d_min; i < p; ++i) {
      if (ncorrelation[k][i] > 0) {
        t[im] = i * pow((double)m, k);
        f[im] = correlation[k][i] / ncorrelation[k][i] - aux;
        ++im;
      }
    }
  }

  npcorr = im;
}
