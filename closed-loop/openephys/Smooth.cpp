#include <math.h>
#include "Smooth.h"

Smooth::Smooth(): m_tau(100), m_srate(1) {
    calc_k();
}

void Smooth::clear() {
    m_smoothed = 0;
}

void Smooth::setTau(float tau) {
    m_tau = tau;
    calc_k();
}

void Smooth::setSampleRate(float srate) {
    m_srate = srate;
    calc_k();
}

void Smooth::calc_k() {
    m_k = exp(-2 * M_PI / ((m_tau / 1000) * m_srate));
    m_smoothed = 0;
}

void Smooth::addSample(double sample) {
    m_smoothed = (m_k * m_smoothed) + ((1 - m_k) * sample);
}

double Smooth::getSmoothed() {
    return m_smoothed;
}
