#ifndef SMOOTH_H
#define SMOOTH_H

#define _USE_MATH_DEFINES

class Smooth {
public:
    Smooth();

    void setTau(float tau);
    void setSampleRate(float srate);
    void addSample(double sample);
    double getSmoothed();
    void clear();

private:
    void calc_k();

    double m_smoothed;

    float m_tau;
    float m_srate;
    float m_k;
};


#endif /* SMOOTH_H */
