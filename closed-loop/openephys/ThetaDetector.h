// Theta ratio calculator
//
// Jinhee Baek


#ifndef THETA_DETECTOR_H_INCLUDED
#define THETA_DETECTOR_H_INCLUDED

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

#include <ProcessorHeaders.h>
#include <DspLib.h>  // Filtering
#include <random>
#include "RMS.h"
#include "ThetaMonitor.h"
#include "ThetaEstimator.h"
#include "Smooth.h"

using BandpassFilter = Dsp::SimpleFilter
            <Dsp::Butterworth::BandPass // filter type
            <2>,                        // order
            1,                          // number of channels
            Dsp::DirectFormII>;         // realization

enum {
	pInputChan,                 // 0-based
	pThetaEventChan,            // 0-based
	pLightEventChan,            // 0-based
	pLightGateChan,             // 0-based
	pThetaLow,
	pThetaHigh,
	pRefLow,
	pRefHigh,
    pThreshold,
    pDelayEnabled
};


class ThetaDetector : public GenericProcessor {
    friend class ThetaDetectorEditor;

public:
    ThetaDetector();
    ~ThetaDetector();

    bool hasEditor() const { return true; }
    AudioProcessorEditor* createEditor() override;

	void updateSettings() override;
	void setFilterParams();

    void process(AudioSampleBuffer& continuousBuffer) override;

    void setParameter(int parameterIndex, float newValue) override;

	bool enable() override;
    bool disable() override;

    void createEventChannels() override;
    void handleEvent(const EventChannel* eventInfo, const MidiMessage& event,
                     int samplePosition = 0) override;

private:
    void triggerThetaCrossing(int tsInBuffer, uint8 direction);
    void triggerLight(int tsInBuffer, int val);
    // void triggerLightOff(int tsInBuffer);
    void setOutputDataArray();

    // input and output channel numbers
    int m_inputChan;
    int m_thetaEventChan;
    int m_lightEventChan;
    int m_lightGateChan;

    // freq ranges
	float m_thetaLow;
	float m_thetaHigh;
	float m_refLow;
	float m_refHigh;

    // temporary buffers and corresponding filters
	AudioSampleBuffer m_thetaBuffer;
	AudioSampleBuffer m_refBuffer;
    BandpassFilter m_thetaFilter;
    BandpassFilter m_refFilter;

    RMS m_rms_ref;
    ThetaEstimator m_theta_estimator;
    Smooth m_smooth_ref;
    Smooth m_smooth_ratio;

    // holded sample during artifact
    float m_sampHoldDur;
    double m_sampHoldTheta;
    double m_sampHoldRef;

    ThetaMonitor m_thetaMonitor;
    float m_threshold;
    float m_tIgnoreThetaAfterLight;

    // light control
    bool m_isLightOn;
    bool m_isDelayEnabled;
    // -1 if light is disabled, >0 if light is enabled from stored time
    juce::int64 m_tsLightSessionEnabled; // assumes that it does not change more than onece within a buffer
    juce::int64 m_tsDelayLight;
    juce::int64 m_tsLightOff;
    juce::int64 m_tsLightOn;

    // External pulse generator should be used for precied duration (due to buffer operation)
    // Fix this value to 0.1 s
    float m_lightDur;

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 m_gen; // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> m_delay;

    // event generation
    EventChannel* m_thetaEventChannelPtr;
    EventChannel* m_lightEventChannelPtr;
    MetaDataDescriptorArray m_thetaEventMetaDataDescriptors;
    MetaDataDescriptorArray m_lightEventMetaDataDescriptors;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ThetaDetector);
};

#endif
