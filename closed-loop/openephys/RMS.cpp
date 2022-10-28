#include "RMS.h"

RMS::RMS() {
	setSize(1000);
}

void RMS::clear() {
    setSize(m_buffer.size());
}

void RMS::setSize(int nsamp) {
	m_buffer.clear();
    m_buffer.insertMultiple(0, 0, nsamp);
	m_index = 0;
	m_ss = 0;
}

void RMS::addSample(double sample) {
    double sq = sample * sample;

    m_ss -= m_buffer[m_index];
	m_ss += sq;

	m_buffer.set(m_index, sq);

	m_index += 1;
	m_index %= m_buffer.size();
}

double RMS::getRMS() {
    return sqrt(m_ss / m_buffer.size());
}
