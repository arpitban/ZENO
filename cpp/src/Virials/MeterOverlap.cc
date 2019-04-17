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

#include "MeterOverlap.h"
#include "alloc2D.h"

///
///

template <class T, class RandomNumberGenerator>
MeterOverlap<T, RandomNumberGenerator>::
MeterOverlap(ClusterSum<T, RandomNumberGenerator> & clusterSumPrimary, ClusterSum<T, RandomNumberGenerator> & clusterSumPerturb, double alphaCenter, double alphaSpan, int numAlpha) : clusterSumPrimary(clusterSumPrimary), clusterSumPerturb(clusterSumPerturb), numAlpha(numAlpha){
    setAlpha(alphaCenter, alphaSpan, numAlpha);
}

template <class T, class RandomNumberGenerator>
MeterOverlap<T, RandomNumberGenerator>::
  ~MeterOverlap() {
    delete[] data;
    delete[] alpha;
}

/// Defines particles using assembly of spheres with mutable center and orientation.
///
template <class T, class RandomNumberGenerator>
void
MeterOverlap<T, RandomNumberGenerator>::
setAlpha(double alphaCenter, double alphaSpan, int nAlpha) {
    numAlpha = nAlpha;
    if (numAlpha > 1) {
        numData = numAlpha;
        if (alphaSpan == 0) {
            std::cerr << "If # of alpha > 1, then alpha span can't be 0" << std::endl;
            exit(1);
        }
    }
    else {
        numData = 2;
        if (alphaSpan != 0) {
            std::cerr << "If # of alpha is 1, then alpha span must be 0" << std::endl;
            exit(1);
        }
    }
    data = new double[numData];
    alpha = new double[numData];

    if (numAlpha == 1) {
        alpha[0] = alphaCenter;
        return;
    }
    for (int i = 0; i < numAlpha; ++i) {
        alpha[i] = alphaCenter * exp((i - (numAlpha - 1.0) / 2) / (numAlpha - 1.0) * alphaSpan);
    }
}

template <class T, class RandomNumberGenerator>
int
MeterOverlap<T, RandomNumberGenerator>::
getNumAlpha() {
 return numAlpha;
}

template <class T, class RandomNumberGenerator>
const double *
MeterOverlap<T, RandomNumberGenerator>::
getAlpha(){
    return alpha;
}

template <class T, class RandomNumberGenerator>
void
MeterOverlap<T, RandomNumberGenerator>::
collectData() {
    const double primaryValue = clusterSumPrimary.value();
    double pi = fabs(primaryValue);
    if (pi == 0 || pi == std::numeric_limits<double>::infinity() || std::isnan(pi)) {
        std::cerr << "pi is" << pi << std::endl;
        exit(1);
    }
    double perturbValue = fabs(clusterSumPerturb.value());
    if (numAlpha == 1) {
        data[0] = primaryValue / pi;
        data[1] = perturbValue / (perturbValue + alpha[0] * pi);
    } else {
        for (int i = 0; i < numAlpha; ++i) {
            // gamma_OS = pi1 pi0 / (pi1 + alpha pi0)
            // 0: gamma_OS/pi0 = pi1 / (pi1 + alpha pi0)
            // 1: gamma_OS/pi1 = pi0 / (pi1 + alpha pi0)
            //                 = (1/alpha) pi0 / (pi1 + (1/alpha) pi0)
            // for 1 case, we use negative alphaSpan (alpha => 1/alpha)
            //    and effectively compute: alpha gammaOS/pi1
            // <0>/<1> = (1/alpha) <gammaOS/pi0>0 / <gammaOS/pi1>1
            //         ~= 1 (when alpha is optimal)
            data[i] = perturbValue / (perturbValue + alpha[i] * pi);
        }
    }
    for (int i = 0; i < numData; ++i) {
        mostRecent[i] = data[i];
        currentBlockSum[i] += data[i];
    }
    if (--blockCountdown == 0) {
        if (doCovariance) {
            double blockSizeSq = ((double) blockSize) * ((double) blockSize);
            for (int i = 0; i < numData; ++i) {
                for (int j = 0; j <= i; ++j) {
                    double ijx = currentBlockSum[i] * currentBlockSum[j] / blockSizeSq;
                    blockCovSum[i][j] += ijx;
                }
            }
        }
        for (int i = 0; i < numData; ++i) {
            blockSum[i] += currentBlockSum[i];
            currentBlockSum[i] /= blockSize;
            if (maxBlockCount > 0) {
                if (blockCount > 0) correlationSum[i] += blockSums[i][blockCount - 1] * currentBlockSum[i];
                blockSums[i][blockCount] = currentBlockSum[i];
            } else {
                if (blockCount > 0) correlationSum[i] += prevBlockSum[i] * currentBlockSum[i];
                else firstBlockSum[i] = currentBlockSum[i];
                prevBlockSum[i] = currentBlockSum[i];
            }
            blockSum2[i] += currentBlockSum[i] * currentBlockSum[i];
            currentBlockSum[i] = 0;
        }
        blockCount++;
        blockCountdown = blockSize;
    }
}

template <class T, class RandomNumberGenerator>
void
MeterOverlap<T, RandomNumberGenerator>::
setBlockSize(long bs) {
    defaultBlockSize = blockSize = bs;
    maxBlockCount = -1;
    reset();
}

template <class T, class RandomNumberGenerator>
void
MeterOverlap<T, RandomNumberGenerator>::
setMaxBlockCount(long mBC) {
    if (mBC < 4) {
        std::cerr << "Max block count needs to be at least 3!"<< std::endl;
        exit(1);
    }
    blockSums = (double ** ) realloc2D ((void ** ) blockSums, numData, maxBlockCount, sizeof(double));
}

template <class T, class RandomNumberGenerator>
void
MeterOverlap<T, RandomNumberGenerator>::
reset() {
    blockCount = 0;
    blockSize = defaultBlockSize;
    blockCountdown = blockSize;
    if (numData == 0) {
        // we can't allocate 0-size arrays, so just leave them as nullptr
        // at some point numData will be positive, reset will be called again
        return;
    }
    // realloc our arrays so that we can adjust if n changes
    mostRecent = (double *) realloc(mostRecent, numData * sizeof(double));
    currentBlockSum = (double *) realloc(currentBlockSum, numData * sizeof(double));
    blockSum = (double *) realloc(blockSum, numData * sizeof(double));
    blockSum2 = (double *) realloc(blockSum2, numData * sizeof(double));
    correlationSum = (double *) realloc(correlationSum, numData * sizeof(double));
    for (int i = 0; i < numData; ++i) {
        currentBlockSum[i] = blockSum[i] = blockSum2[i] = correlationSum[i] = 0;
    }
    if (maxBlockCount > 0) {
        if (maxBlockCount % 2 == 1 || maxBlockCount < 4) {
            std::cerr << "Not nice!  Give me a max block count that's even and >= 4!" << std::endl;
            exit(1);
        }
        blockSums = (double **) realloc2D((void **) blockSums, numData, maxBlockCount, sizeof(double));
    } else {
        prevBlockSum = (double *) realloc(prevBlockSum, numData * sizeof(double));
        firstBlockSum = (double *) realloc(firstBlockSum, numData * sizeof(double));
    }
    stats = (double **) realloc2D((void **) stats, numData, 4, sizeof(double));
    if (doCovariance) {
        blockCovSum = (double **) realloc2D((void **) blockCovSum, numData, numData, sizeof(double));
        for (int i = 0; i < numData; ++i) {
            for (int j = 0; j <= i; ++j) {
                blockCovSum[i][j] = 0;
            }
        }
        blockCovariance = (double **) realloc2D((void **) blockCovariance, numData, numData, sizeof(double));
    }
    ratioStats = (double ** ) realloc2D ((void ** ) ratioStats, numData, 3, sizeof(double));
    if (doCovariance) ratioCovariance = (double ** ) realloc2D ((void ** ) ratioCovariance, numData, numData, sizeof(double));
}

template <class T, class RandomNumberGenerator>
double **
MeterOverlap<T, RandomNumberGenerator>::
getStatistics() {
    if (blockCount == 0) {
        for (int i = 0; i < numData; ++i) {
            stats[i][AVG_CUR] = mostRecent[i];
            stats[i][1] = stats[i][2] = stats[i][3] = NAN;
        }
        return stats;
    }
    for (int i = 0; i < numData; ++i) {
        stats[i][AVG_CUR] = mostRecent[i];
        stats[i][AVG_AVG] = blockSum[i] / (blockSize * blockCount);
        if (blockCount == 1) {
            for (int i = 0; i < numData; ++i) {
                stats[i][AVG_ERR] = stats[i][AVG_ACOR] = NAN;
            }
            continue;
        }
        stats[i][AVG_ERR] = blockSum2[i] / blockCount - stats[i][AVG_AVG] * stats[i][AVG_AVG];
        if (stats[i][AVG_ERR]<0) stats[i][AVG_ERR] = 0;
        if (stats[i][AVG_ERR] == 0) {
            stats[i][AVG_ACOR] = 0;
        }
        else {
            double bc;
            if (maxBlockCount>0) {
                bc = (((2 * blockSum[i] / blockSize - blockSums[i][0] - blockSums[i][blockCount - 1]) * stats[i][AVG_AVG] - correlationSum[i]) / (1 - blockCount) + stats[i][AVG_AVG] * stats[i][AVG_AVG]) / stats[i][AVG_ERR];
            }
            else {
                bc = (((2 * blockSum[i] / blockSize - firstBlockSum[i] - prevBlockSum[i]) * stats[i][AVG_AVG] - correlationSum[i]) / (1 - blockCount) + stats[i][AVG_AVG] * stats[i][AVG_AVG]) / stats[i][AVG_ERR];
            }
            stats[i][AVG_ACOR] = (std::isnan(bc) || bc <= -1 || bc >= 1) ? 0 : bc;
        }
        stats[i][AVG_ERR] = sqrt(stats[i][AVG_ERR] / (blockCount - 1));
    }
    return stats;
}

template <class T, class RandomNumberGenerator>
double **
MeterOverlap<T, RandomNumberGenerator>::
getBlockCovariance() {
    double totSamples = blockSize*blockCount;
    double totSq = totSamples*totSamples;
    for (int i = 0; i < numData; ++i) {
        for (int j = 0; j <= i; ++j) {
            blockCovariance[i][j] = blockCovariance[j][i] = blockCovSum[i][j] / blockCount - blockSum[i] * blockSum[j] / totSq;
        }
    }

    return blockCovariance;
}

template <class T, class RandomNumberGenerator>
double **
MeterOverlap<T, RandomNumberGenerator>::
getBlockCorrelation() {
    getBlockCovariance();
    for (int i = 0; i < numData; ++i) {
        for (int j = 0; j < i; ++j) {
            double d = blockCovariance[i][i] * blockCovariance[j][j];
            double c = d <= 0 ? 0 : blockCovariance[i][j] / sqrt(d);
            blockCovariance[i][j] = blockCovariance[j][i] = c;
        }
    }
    for (int i = 0; i < numData; ++i) {
        blockCovariance[i][i] = 1;
    }
    return blockCovariance;
}

template <class T, class RandomNumberGenerator>
double **
MeterOverlap<T, RandomNumberGenerator>::
getRatioStatistics() {
    if (blockCount == 0) {
        for (int i = 0; i < numData; ++i) {
            ratioStats[i][AVG_CUR] = NAN;
            ratioStats[i][AVG_AVG] = ratioStats[i][AVG_ERR] = NAN;
        }
        return ratioStats;
    }
    getStatistics();
    getBlockCovariance();
    for (int i = 0; i < numData; ++i) {
        ratioStats[i][AVG_CUR] = mostRecent[i] / mostRecent[numData - 1];
        ratioStats[i][AVG_AVG] = blockSum[i] / blockSum[numData - 1];
        if (blockCount == 1) {
            for (int i = 0; i < numData; ++i) {
                ratioStats[i][AVG_ERR] = NAN;
            }
            continue;
        }
        double d = blockCovariance[i][i] * blockCovariance[numData - 1][numData - 1];
        double icor = d <= 0 ? 0 : blockCovariance[i][numData - 1] / sqrt(d);
        ratioStats[i][AVG_ERR] = ratioErr(stats[i][AVG_AVG], stats[i][AVG_ERR], stats[numData - 1][AVG_AVG], stats[numData - 1][AVG_ERR], icor);
    }
    return ratioStats;
}

template <class T, class RandomNumberGenerator>
double **
MeterOverlap<T, RandomNumberGenerator>::
getRatioCovariance() {
    if (blockCount<2) {
        for (int i = 0; i < numData; ++i) {
            for (int j = 0; j < numData; ++j) {
                ratioCovariance[i][j] = NAN;
            }
        }
        return ratioCovariance;
    }
    getStatistics();
    getBlockCovariance();
    double vd = stats[numData - 1][AVG_AVG];
    double ed = stats[numData - 1][AVG_ERR];
    for (int i = 0; i < numData; ++i) {
        double vi = stats[i][AVG_AVG];
        double ei = stats[i][AVG_ERR];
        double x = blockCovariance[i][i] * blockCovariance[numData - 1][numData - 1];
        if (x <= 0) {
            for (int j = 0; j <= i; ++j) ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
            continue;
        }
        double cid = blockCovariance[i][numData - 1] / sqrt(x);
        for (int j = 0; j <= i; ++j) {
            double vj = stats[j][AVG_AVG];
            double ej = stats[j][AVG_ERR];
            if (blockCovariance[j][j] <= 0) {
                ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
                continue;
            }
            double cjd = blockCovariance[j][numData - 1] / sqrt(blockCovariance[j][j] * blockCovariance[numData - 1][numData - 1]);
            double cij = blockCovariance[i][j] / sqrt(blockCovariance[i][i]*blockCovariance[j][j]);
            ratioCovariance[j][i] = ratioCovariance[i][j] = ratioCov(vi, vj, vd, ei, ej, ed, cij, cid, cjd);
        }
    }
    return ratioCovariance;
}

template <class T, class RandomNumberGenerator>
double **
MeterOverlap<T, RandomNumberGenerator>::
getRatioCorrelation() {
    if (blockCount<2) {
        for (int i = 0; i < numData; ++i) {
            for (int j = 0; j < numData; ++j) {
                ratioCovariance[i][j] = NAN;
            }
        }
        return ratioCovariance;
    }
    getStatistics();
    getBlockCovariance();
    double vd = stats[numData - 1][AVG_AVG];
    double ed = stats[numData - 1][AVG_ERR];
    for (int i = 0; i < numData; ++i) {
        double vi = stats[i][AVG_AVG];
        double ei = stats[i][AVG_ERR];
        double x = blockCovariance[i][i] * blockCovariance[numData - 1][numData - 1];
        if (x <= 0) {
            for (int j = 0; j <= i; ++j) ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
            continue;
        }
        double cid = blockCovariance[i][numData - 1] / sqrt(x);
        for (int j = 0; j <= i; ++j) {
            double vj = stats[j][AVG_AVG];
            double ej = stats[j][AVG_ERR];
            if (blockCovariance[j][j] <= 0) {
                ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
                continue;
            }
            double cjd = blockCovariance[j][numData - 1] / sqrt(blockCovariance[j][j] * blockCovariance[numData - 1][numData - 1]);
            double cij = blockCovariance[i][j] / sqrt(blockCovariance[i][i]*blockCovariance[j][j]);
            ratioCovariance[j][i] = ratioCovariance[i][j] = ratioCor(vi, vj, vd, ei, ej, ed, cij, cid, cjd);
        }
    }
    return ratioCovariance;
}