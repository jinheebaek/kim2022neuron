#include <cmath>
#include "ThetaMonitor.h"

ThetaMonitor::ThetaMonitor()
    : m_threshold(0.5), m_sampleRate(1000), m_twindow(0.2),
      m_isUp(false), m_isPrevUp(false),
      m_index(0), m_count(0) {
}

void ThetaMonitor::clear() {
    resetWindow();
}

bool ThetaMonitor::isUp() { return m_isUp; }

bool ThetaMonitor::isCrossingUp() { return (!m_isPrevUp && m_isUp); }

bool ThetaMonitor::isCrossingDown() { return (m_isPrevUp && !m_isUp); }

void ThetaMonitor::setSampleRate(float sampleRate) {
    m_sampleRate = sampleRate;
    resetWindow();
}

void ThetaMonitor::setTwindow(float twindow) {
    m_twindow = twindow;
    resetWindow();
}

void ThetaMonitor::setThreshold(float threshold) {
    m_threshold = threshold;
}

void ThetaMonitor::resetWindow() {
    int size = static_cast<int>(std::floor(m_sampleRate * m_twindow));
    m_index = 0;
    m_count = 0;
	m_data.clear();
	m_data.insertMultiple(0, false, size);
    m_isUp = false;
    m_isPrevUp = false;
}

void ThetaMonitor::addSample(double sample) {
    bool result = sample > m_threshold;

    m_count -= m_data[m_index];
	m_count += result;

    m_data.set(m_index, result);

    m_index += 1;
    m_index %= m_data.size();

    m_isPrevUp = m_isUp;
    if (!m_isUp)
        m_isUp = (m_count == m_data.size());
    else
        m_isUp = (m_count > 0);
}
