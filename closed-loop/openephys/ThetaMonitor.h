#ifndef THETAMONITOR_H
#define THETAMONITOR_H

#include <ProcessorHeaders.h>

class ThetaMonitor {
public:
    ThetaMonitor();
    bool isUp();
    bool isCrossingUp();
    bool isCrossingDown();

    void addSample(double sample);
    void clear();

    void setSampleRate(float sampleRate);
    void setTwindow(float twindow);
    void setThreshold(float threshold);

private:
    void resetWindow();

    int m_index;
    juce::Array<bool> m_data;
    double m_count;

    bool m_isPrevUp;
    bool m_isUp;

    float m_threshold;
    float m_sampleRate;
    float m_twindow;
};

#endif /* THETAMONITOR_H */
