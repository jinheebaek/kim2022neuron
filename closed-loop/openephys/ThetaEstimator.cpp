#include <assert.h>
#include <algorithm>
#include "ThetaEstimator.h"
#include <iostream>

ThetaEstimator::ThetaEstimator():
    m_thetaLow(4), m_thetaHigh(8), m_bw_sub(1) {
    reset();
}

void ThetaEstimator::clear() {
    for (auto& smooth : m_smooth_subbands)
        smooth.clear();
}

void ThetaEstimator::setThetaLow(float thetaLow) {
    m_thetaLow = thetaLow;
    reset();
}

void ThetaEstimator::setThetaHigh(float thetaHigh) {
    m_thetaHigh = thetaHigh;
    reset();
}

void ThetaEstimator::setSubBandWidth(float bw) {
    m_bw_sub = bw;
    reset();
}

void ThetaEstimator::setSampleRate(float srate) {
    SDFT::setSampleRate(srate);
}

void ThetaEstimator::reset() {
    float bw_theta = m_thetaHigh - m_thetaLow;

    assert(bw_theta >= m_bw_sub);

    m_nsubband = std::ceil((bw_theta - m_bw_sub) / (m_bw_sub / 2)) + 1;
    m_fstep = (bw_theta - m_bw_sub) / (m_nsubband - 1);

    m_smooth_subbands.clear();
    m_smooth_subbands.reserve(m_nsubband);
    for (int i = 0; i < m_nsubband; i++) {
        m_smooth_subbands.push_back(Smooth());
        m_smooth_subbands[i].setTau(250);
    }
}

void ThetaEstimator::addSample(double sample) {
    SDFT::addSample(sample);

    auto it = m_smooth_subbands.begin();
    for (int i = 0; i < m_nsubband; i++) {
        float fmin = m_thetaLow + (i * m_fstep);
        float fmax = fmin + m_bw_sub;

        it->addSample(getBandPower(fmin, fmax));
    }
}

double ThetaEstimator::getThetaPower() {
    std::vector<double> subpower;
    subpower.reserve(m_nsubband);
    for (auto& smooth : m_smooth_subbands)
        subpower.push_back(smooth.getSmoothed());

    return *std::max_element(subpower.begin(), subpower.end());
}
