// ================================================================
//
// This software was developed by employees of the National Institute of
// Standards and Technology (NIST), an agency of the Federal Government.
// Pursuant to title 17 United States Code Section 105, works of NIST employees
// are not subject to copyright protection in the United States and are
// considered to be in the public domain. Permission to freely use, copy,
// modify, and distribute this software and its documentation without fee is
// hereby granted, provided that this notice and disclaimer of warranty appears
// in all copies.
//
// THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
// EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
// WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
// FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
// THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
// EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
// DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
// RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
// BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
// SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
// SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
// SERVICES PROVIDED HEREUNDER.
//
// Distributions of NIST software should also include copyright and licensing
// statements of any third-party software that are legally bundled with the
// code in compliance with the conditions of those licenses.
// 
// ================================================================

// ================================================================
// 
// Authors:
// Created:
//
// ================================================================

#ifndef METER_OVERLAP_H
#define METER_OVERLAP_H

#include "ClusterSum.h"

///
///
template <class T,
        class RandomNumberGenerator>
class ClusterSum;

template <class T, class RandomNumberGenerator>
class MeterOverlap {
 public:
    MeterOverlap(ClusterSum<T, RandomNumberGenerator> & clusterSumPrimary, ClusterSum<T, RandomNumberGenerator> & clusterSumPerturb, double alphaCenter, double alphaSpan, int numAlpha);
    ~MeterOverlap();
    void setAlpha(double alphaCenter, double alphaSpan, int nAlpha);
    int getNumAlpha();
    const double * getAlpha();
    void collectData();
    double ** getStatistics();
    double ** getBlockCovariance();
    double ** getBlockCorrelation();
    void setBlockSize(long blockSize);
    void setMaxBlockCount(long maxBlockCount);
    long getBlockSize() {return blockSize;}
    long getBlockCount() {return blockCount;}
    void setNumData(int newNumData);
    int getNumData();
    virtual void reset();
    double ** getRatioStatistics();
    double ** getRatioCovariance();
    double ** getRatioCorrelation();
    const static int AVG_CUR = 0, AVG_AVG = 1, AVG_ERR = 2, AVG_ACOR = 3;

    static double ratioErr(double nAvg, double nErr, double dAvg, double dErr, double cor) {
        if (nAvg == 0 && nErr == 0) return 0;
        double ratio = nAvg / dAvg;
        if (nAvg == 0) {
            return sqrt((nErr * nErr) / (dAvg * dAvg));
        }
        return sqrt((nErr * nErr / (nAvg * nAvg) + dErr * dErr / (dAvg * dAvg) - 2 * cor * nErr * dErr / (nAvg * dAvg)) * ratio * ratio);
    }
    /**
 * Compute covariance of i/d and j/d
 *
 * vi: value of i
 * vj: value of j
 * vd: value of d
 * ei: error in the i numerator
 * ej: error in the j numerator
 * ed: error in the denominator d
 * cij: correlation between i and j
 * cid: correlation between i and d
 * cjd: correlation between j and d
 */
    static double ratioCov(double vi, double vj, double vd, double ei, double ej, double ed, double cij, double cid, double cjd) {
        double eid = ratioErr(vi, ei, vd, ed, cid);
        double ejd = ratioErr(vj, ej, vd, ed, cjd);
        if (eid == 0 || ejd == 0) return 0;
        return (vi / vd) * (vj / vd) * ((ei / vi) * (ed / vd) + (ei / vi) * (ej / vj) * cij - (ei / vi) * (ed / vd) * cid - (ej / vj) * (ed / vj) * cjd);
    }

    /**
     * Compute correlation of i/d and j/d
     *
     * vi: value of i
     * vj: value of j
     * vd: value of d
     * ei: error in the i numerator
     * ej: error in the j numerator
     * ed: error in the denominator d
     * cij: correlation between i and j
     * cid: correlation between i and d
     * cjd: correlation between j and d
     */
    static double ratioCor(double vi, double vj, double vd, double ei, double ej, double ed, double cij, double cid, double cjd) {
        double eid = ratioErr(vi, ei, vd, ed, cid);
        double ejd = ratioErr(vj, ej, vd, ed, cjd);
        if (eid == 0 || ejd == 0) return 0;
        return ((vi / vd) / eid) * ((vj / vd) / ejd) * ((ed / vd) * (ed / vd) + (ei / vi) * (ej / vj) * cij - (ei / vi) * (ed / vd) * cid - (ej / vj) * (ed / vd) * cjd);
    }
private:
    ClusterSum<T, RandomNumberGenerator> & clusterSumPrimary;
    ClusterSum<T, RandomNumberGenerator> & clusterSumPerturb;
    double * data;
    double * alpha;
    const int numAlpha;
    int numData;
    long defaultBlockSize, blockSize, blockCount, maxBlockCount;
    long blockCountdown;
    double * mostRecent;
    double * currentBlockSum, * blockSum, * blockSum2, * correlationSum;
    double * prevBlockSum, * firstBlockSum;
    double ** stats;
    double ** blockSums;
    double ** blockCovariance;
    double ** blockCovSum;
    const bool doCovariance;
    double ** ratioStats;
    double ** ratioCovariance;
};
#endif

