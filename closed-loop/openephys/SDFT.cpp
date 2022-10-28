#include <math.h>
#include "SDFT.h"
#include <iostream>
#include <algorithm>

SDFT::SDFT(): m_srate(1), m_max_freq(100) {
    setNfft(4096);
}

void SDFT::clear() {
    m_index = 0;
    m_sample.clear();
    m_fftout.clear();
}

void SDFT::setNfft(int nfft) {
    m_nfft = nfft;
    m_index = 0;

    m_downsample_factor = 10;
    m_downsample_counter = 0;

    m_sample.resize(nfft);
    m_sample.assign(0, nfft);
    m_fftout.resize(nfft);
    m_fftout.assign(0, nfft);

    m_coefs.resize(nfft);
    for (int i = 0; i < nfft; ++i) {
        double a = 2.0 * M_PI * i  / static_cast<double>(nfft);
        m_coefs[i] = complex(cos(a)/* / N */, sin(a) /* / N */);
    }

    setSampleRate(m_srate);
}

void SDFT::setSampleRate(float srate) {
    m_srate = srate;

    m_max_idx = std::floor(m_max_freq / (m_srate / m_nfft));
    m_max_idx = m_nfft < m_max_idx ? m_nfft : m_max_idx;
}

void SDFT::addSample(double sample) {
    m_downsample_counter++;
    m_downsample_counter %= m_downsample_factor;

    // if (!m_downsample_counter == 0)
    //     return;

    complex delta = sample - m_sample[m_index];
    m_sample[m_index] = sample;

    // for (int i = 0; i < m_nfft; ++i) {
    for (int i = 0; i < m_max_idx; ++i) {
        m_fftout[i] = (m_fftout[i] + delta) * m_coefs[i];
    }

    m_index++;
    m_index %= m_nfft;
}

double SDFT::getBandPower(float fmin, float fmax) {
    double power = 0;
    int imin = std::floor(fmin / (m_srate / m_nfft)) + 1;
    int imax = std::floor(fmax / (m_srate / m_nfft));

    assert(imin <= m_max_idx);
    assert(imax <= m_max_idx);

    for (int i = imin; i <= imax; i++) {
        complex v = m_fftout[i];
        power += v.real() * v.real() + v.imag() * v.imag();
    }
    return power / m_nfft;
}
