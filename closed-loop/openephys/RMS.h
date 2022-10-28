#ifndef RMS_H
#define RMS_H

#include <ProcessorHeaders.h>

class RMS {
public:
    RMS();
	void setSize(int nsamp);
    void addSample(double sample);
    double getRMS();
    void clear();

private:
    juce::Array<double> m_buffer;     // juce array
	int m_index;
	double m_ss;
};


#endif /* RMS_H */
