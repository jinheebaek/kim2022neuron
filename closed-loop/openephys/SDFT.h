#ifndef SDFT_H
#define SDFT_H

#define _USE_MATH_DEFINES

#include <vector>
#include <complex>
#include "Smooth.h"

typedef std::complex<double> complex;

class SDFT {
public:
    SDFT();

    void setNfft(int nfft);
    void setSampleRate(float srate);
    void clear();

    virtual void addSample(double sample);
    double getBandPower(float fmin, float fmax);

private:
    int m_nfft;
    int m_index;
    float m_srate;

    float m_max_freq;
    float m_max_idx;

    int m_downsample_factor;
    int m_downsample_counter;

    // SDFT
    std::vector<complex> m_sample;
    std::vector<complex> m_coefs;
    std::vector<complex> m_fftout;
};


#endif /* SDFT_H */
