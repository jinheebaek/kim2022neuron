#ifndef ThetaEstimator_H
#define ThetaEstimator_H

#include "SDFT.h"
#include "Smooth.h"

class ThetaEstimator: SDFT {
public:
    ThetaEstimator();

    void addSample(double sample) override;
    double getThetaPower();
    void clear();

    void setThetaLow(float thetaLow);
    void setThetaHigh(float thetaHigh);
    void setSubBandWidth(float bw);
    void setSampleRate(float srate);

private:
    void reset();

    float m_thetaLow;
    float m_thetaHigh;
    float m_bw_sub;

    std::vector<Smooth> m_smooth_subbands;
    int m_nsubband;
    float m_fstep;
};


#endif /* ThetaEstimator_H */
