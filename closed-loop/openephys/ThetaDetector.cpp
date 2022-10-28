/*
------------------------------------------------------------------

This file is part of a plugin for the Open Ephys GUI
Copyright (C) 2017 Translational NeuroEngineering Laboratory, MGH

------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
/*
This is based heavily on the crossing detector third party plugin
First, simply filtering the data before thresholding
*/

#include <math.h>
#include <random>
#include "ThetaDetector.h"
#include "ThetaDetectorEditor.h"

ThetaDetector::ThetaDetector()
    : GenericProcessor("Theta Detector"), m_inputChan(0), m_thetaEventChan(0),
      m_lightEventChan(0), m_lightGateChan(2),
      m_thetaLow(4.0f), m_thetaHigh(8.0f), m_refLow(1.0f), m_refHigh(20.0f),
      m_sampHoldTheta(0.0f), m_sampHoldRef(0.0f),
      m_threshold(0.5), m_isLightOn(false), m_tIgnoreThetaAfterLight(1),
      m_isDelayEnabled(false), m_lightDur(0.1f), m_sampHoldDur(0.05f),
      m_tsLightSessionEnabled(-1.0f), m_tsDelayLight(-1.0f),
      m_gen(rd()), m_delay(0.5, 1.0)      {
    setProcessorType(PROCESSOR_TYPE_FILTER);

    // make the event-related metadata descriptors
    // theta event
    m_thetaEventMetaDataDescriptors.add(new MetaDataDescriptor(MetaDataDescriptor::INT64, 1, "Crossing Point",
        "Time when theta threshold was crossed", "crossing.point"));
    m_thetaEventMetaDataDescriptors.add(new MetaDataDescriptor(MetaDataDescriptor::FLOAT, 1, "Threshold",
        "Theta ratio threshold", "crossing.threshold"));
    m_thetaEventMetaDataDescriptors.add(new MetaDataDescriptor(MetaDataDescriptor::UINT8, 1, "Direction",
        "Direction of crossing: 1 = rising, 0 = falling", "crossing.direction"));

    // light event
    m_lightEventMetaDataDescriptors.add(new MetaDataDescriptor(MetaDataDescriptor::INT64, 1, "Light Event Point",
        "Time when light was triggered", "light.point"));
    // OFF event does not represent actual light duration
    // here, duration indicate duration of light triggering stimulation
    // currently light trigger stim dur was fixed to 0.1
    m_lightEventMetaDataDescriptors.add(new MetaDataDescriptor(MetaDataDescriptor::UINT8, 1, "Direction",
        "Direction of stim: 1 = on, 0 = off", "light.direction"));
}

ThetaDetector::~ThetaDetector() {
}

AudioProcessorEditor* ThetaDetector::createEditor() {
    editor = new ThetaDetectorEditor(this);
    return editor;
}

void ThetaDetector::setOutputDataArray() {
    if (dataChannelArray.size() > 0) {
        OwnedArray<DataChannel> upstreamDataChannelArray;
        upstreamDataChannelArray.swapWith(dataChannelArray);

        float srate = upstreamDataChannelArray[0]->getSampleRate();
        dataChannelArray.clear();

        DataChannel* rawch = new DataChannel(*upstreamDataChannelArray[0]);
        dataChannelArray.add(rawch);

        DataChannel* thetach = new DataChannel(*upstreamDataChannelArray[0]);
        thetach->setName("theta");
        dataChannelArray.add(thetach);

        DataChannel* ratioch = new DataChannel(*upstreamDataChannelArray[0]);
        ratioch->setName("ratio");
        dataChannelArray.add(ratioch);

        settings.numOutputs = 3;
    }
}

void ThetaDetector::updateSettings() {
	if (settings.numInputs > 0) {
        setOutputDataArray();

        // const DataChannel* in = upstreamDataChannelArray[m_inputChan];
        const DataChannel* in = getDataChannel(0);
		float sampleRate = in ? in->getSampleRate() : CoreServices::getGlobalSampleRate();

		m_thetaBuffer = AudioSampleBuffer(1, sampleRate);
		m_refBuffer = AudioSampleBuffer(1, sampleRate);

		setFilterParams();

        m_thetaMonitor.setSampleRate(sampleRate);
        m_thetaMonitor.setThreshold(m_threshold);

        m_theta_estimator.setSampleRate(sampleRate);
        int rmsBuffSize = static_cast<int>(sampleRate * 0.5);
        // m_rms_theta.setSize(rmsBuffSize);
        m_rms_ref.setSize(rmsBuffSize);

        m_smooth_ref.setSampleRate(sampleRate);
        m_smooth_ref.setTau(250);
        m_smooth_ratio.setSampleRate(sampleRate);
        m_smooth_ratio.setTau(250);
	}
}

void ThetaDetector::setFilterParams() {
	if (settings.numInputs == 0)
		return;

	//design 3 filters with similar properties
	// int sampRate = upstreamDataChannelArray[m_inputChan]->getSampleRate();
	int sampRate = getDataChannel(0)->getSampleRate();

    m_thetaFilter.setup(
        2,                              // order
        sampRate,                       // sample rate
        (m_thetaHigh + m_thetaLow) / 2, // center frequency
        m_thetaHigh - m_thetaLow        // bandwidth
        );
    m_refFilter.setup(
        2,                          // order
        sampRate,                   // sample rate
        (m_refHigh + m_refLow) / 2, // center frequency
        m_refHigh - m_refLow        // bandwidth
        );
}

bool ThetaDetector::enable() {
    m_tsLightOn = -m_sampHoldDur;
    m_tsLightOff = -m_sampHoldDur;
    m_tsDelayLight = -1;
    m_isLightOn = false;

    m_thetaMonitor.clear();
    m_theta_estimator.clear();
    // m_rms_ref.clear();
    // m_smooth_ref.clear();
    m_smooth_ratio.clear();

    return true;
}

void ThetaDetector::createEventChannels() {
    // add detection event channel
    // const DataChannel* in = upstreamDataChannelArray[m_inputChan];
    const DataChannel* in = getDataChannel(0);

    if (!in) {
        m_thetaEventChannelPtr = nullptr;
        m_lightEventChannelPtr = nullptr;
        return;
    }

    float sampleRate = in->getSampleRate();

    // metadata storing source data channel
    MetaDataDescriptor sourceChanDesc(MetaDataDescriptor::UINT16, 3, "Source Channel",
        "Index at its source, Source processor ID and Sub Processor index of the channel that triggers this event", "source.channel.identifier.full");
    MetaDataValue sourceChanVal(sourceChanDesc);
    uint16 sourceInfo[3];
    sourceInfo[0] = in->getSourceIndex();
    sourceInfo[1] = in->getSourceNodeID();
    sourceInfo[2] = in->getSubProcessorIdx();
    sourceChanVal.setValue(static_cast<const uint16*>(sourceInfo));

    //// theta events
    // EventChannel* thetachan = new EventChannel(EventChannel::TTL, 8, 1, sampleRate, this);
    EventChannel* thetachan = new EventChannel(EventChannel::TTL, 4, 1, sampleRate, this);
    thetachan->setName("Theta crossing output");
    thetachan->setDescription("Triggers whenever the theta ratio crosses a voltage threshold.");
    thetachan->setIdentifier("theta_crossing.event");
    thetachan->addMetaData(sourceChanDesc, sourceChanVal);

    // event-related metadata!
    for (auto desc : m_thetaEventMetaDataDescriptors) {
        thetachan->addEventMetaData(desc);
    }
    m_thetaEventChannelPtr = eventChannelArray.add(thetachan);

    //// light events
    // EventChannel* lightchan = new EventChannel(EventChannel::TTL, 8, 1, sampleRate, this);
    EventChannel* lightchan = new EventChannel(EventChannel::TTL, 4, 1, sampleRate, this);
    lightchan->setName("Light trigger output");
    lightchan->setDescription("Light stimulation for closed-loop manipulation.");
    lightchan->setIdentifier("light.event");
    lightchan->addMetaData(sourceChanDesc, sourceChanVal);

    // event-related metadata!
    for (auto desc : m_lightEventMetaDataDescriptors) {
        lightchan->addEventMetaData(desc);
    }
    m_lightEventChannelPtr = eventChannelArray.add(lightchan);
}

void ThetaDetector::process(AudioSampleBuffer& continuousBuffer) {
   // for open ephys event handling
    checkForEvents();

    // int nSamples = getNumSamples(m_inputChan);
    int nSamples = continuousBuffer.getNumSamples();
	// int sampleRate = upstreamDataChannelArray[m_inputChan]->getSampleRate();
	int sampleRate = getDataChannel(0)->getSampleRate();

    // save raw signal from selected channel
    continuousBuffer.copyFrom(
        0,
        0,
        continuousBuffer,
        m_inputChan,
        0,
        nSamples);

    // copy to buffer
    m_thetaBuffer.copyFrom(0,
                           0,
                           continuousBuffer,
                           m_inputChan,
                           0,
                           nSamples);
    m_refBuffer.copyFrom(0,
                         0,
                         continuousBuffer,
                         m_inputChan,
                         0,
                         nSamples);

	// filter channel in theta band
	float* ptrTw = m_thetaBuffer.getWritePointer(0);
	m_thetaFilter.process(nSamples, &ptrTw);

	// filter copied channel in ref band
	float* ptrRw = m_refBuffer.getWritePointer(0);
	m_refFilter.process(nSamples, &ptrRw);

    const float* ptrT = m_thetaBuffer.getReadPointer(0);
    const float* ptrR = m_refBuffer.getReadPointer(0);

    // int tsBuffer = getTimestamp(m_inputChan);
    int tsBuffer = getTimestamp(0);
    int ts = tsBuffer;
    double p_theta, rms_ref;
    double ratio;
    int tHoldDur = static_cast<int>(sampleRate * m_sampHoldDur);

    for (int i = 0; i < nSamples; i++) {
        if ((tsBuffer + i - m_tsLightOn) <= tHoldDur) {         // if near light stim
            m_theta_estimator.addSample(m_sampHoldTheta);
            m_rms_ref.addSample(m_sampHoldRef);
        } else {
            m_theta_estimator.addSample(*ptrT);
            m_rms_ref.addSample(*ptrR);
        }

        p_theta = m_theta_estimator.getThetaPower();

        m_smooth_ref.addSample(m_rms_ref.getRMS());
        rms_ref = m_smooth_ref.getSmoothed();

        if (rms_ref != 0)
            m_smooth_ratio.addSample(sqrt(p_theta) / rms_ref);
        ratio = m_smooth_ratio.getSmoothed();

        // output theta
        continuousBuffer.setSample(1, i, *ptrT);

        // output theta ratio
        continuousBuffer.setSample(2, i, ratio);

        // for testing
        m_tsLightSessionEnabled = 0;

        m_thetaMonitor.addSample(ratio);
        if (m_thetaMonitor.isCrossingUp() &&
            ((tsBuffer + i - m_tsLightOn) >= m_tIgnoreThetaAfterLight)) {
            triggerThetaCrossing(i, 1); // 1 : over threshold

            // assumes that lightEnable gate does not change more than onece within a buffer length
            if (!m_isDelayEnabled &&                                 // stim. without delay
                (m_tsLightSessionEnabled >= 0) &&                    // light session on
                ((tsBuffer + i - m_tsLightSessionEnabled) >= 0)) {   // light session on
                triggerLight(i, 1);
            }

            m_tsDelayLight = -1; // cancel prev. delayed light, if any.

        } else if (m_thetaMonitor.isCrossingDown()) {
            triggerThetaCrossing(i, 0); // 0 : under threshold

            // reserve delayed light stim.
            if (m_isDelayEnabled) {
                float delay = m_delay(m_gen);
                m_tsDelayLight = tsBuffer + i + static_cast<int>(delay * sampleRate);
            }

        } else {
            // Triggering delayed light after delay
            if ((m_tsDelayLight > 0) &&                              // delay triggered
                (m_tsDelayLight <= (tsBuffer + i)) &&                // delay has passed
                (!m_thetaMonitor.isUp()) &&                          // theta is off
                (m_tsLightSessionEnabled >= 0) &&                    // light session on
                ((tsBuffer + i - m_tsLightSessionEnabled) >= 0)) {   // light session on
                triggerLight(i, 1);
            }
        }

        // light off if duration has passed
        if (m_isLightOn && (tsBuffer + i >= m_tsLightOff)) {
            triggerLight(i, 0);
        }

        // Sample holding right after light stim
        if ((tsBuffer + i) == m_tsLightOn) {
            m_sampHoldTheta = *ptrT;
            m_sampHoldRef = *ptrR;
        }

        ptrT++;
        ptrR++;
    }
}

void ThetaDetector::handleEvent(const EventChannel* eventInfo, const MidiMessage& event, int samplePosition) {
    std::cout << "even detected" << std::endl;

    if (Event::getEventType(event) == EventChannel::TTL) {
		TTLEventPtr ttl = TTLEvent::deserializeFromMessage(event, eventInfo);

        const int eventId         = ttl->getState() ? 1: 0;
        const int eventChannel    = ttl->getChannel();

        if (eventChannel == m_lightGateChan) {
            if (eventId == 1) {
                m_tsLightSessionEnabled = ttl->getTimestamp();
                std::cout << "light session enabled" << std::endl;
            }
            else {
                m_tsLightSessionEnabled = -1;
                std::cout << "light session disabled" << std::endl;
            }
        }
    }
}

void ThetaDetector::triggerThetaCrossing(int tsInBuffer, uint8 direction) {
    // Construct metadata array
    // The order has to match the order the descriptors are stored in createEventChannels.
    MetaDataValueArray mdArray;
    m_thetaEventChan = 3;

    int mdInd = 0;
    MetaDataValue* crossingPointVal = new MetaDataValue(
        *m_thetaEventMetaDataDescriptors[mdInd++]);
    crossingPointVal->setValue(getTimestamp(0) + tsInBuffer);
    mdArray.add(crossingPointVal);

    MetaDataValue* threshVal = new MetaDataValue(
        *m_thetaEventMetaDataDescriptors[mdInd++]);
    threshVal->setValue(m_threshold);
    mdArray.add(threshVal);

    MetaDataValue* directionVal = new MetaDataValue(
        *m_thetaEventMetaDataDescriptors[mdInd++]);
    directionVal->setValue(static_cast<juce::uint8>(direction));
    mdArray.add(directionVal);

    uint8 ttlData;
    if (direction)
        ttlData = 1 << m_thetaEventChan;
    else
        ttlData = 0;

    TTLEventPtr event = TTLEvent::createTTLEvent(
        m_thetaEventChannelPtr,
        getTimestamp(0) + tsInBuffer,
        &ttlData, sizeof(uint8), mdArray, m_thetaEventChan
        );
    addEvent(m_thetaEventChannelPtr, event, tsInBuffer);
}

void ThetaDetector::triggerLight(int tsInBuffer, int val) {
    if ((m_isLightOn && val) || (!m_isLightOn && !val))
        return;

    int ts = getTimestamp(0) + tsInBuffer;

    // Construct metadata array
    // The order has to match the order the descriptors are stored in createEventChannels.
    MetaDataValueArray mdArray;
    int mdInd = 0;

    MetaDataValue* crossingPointVal = new MetaDataValue(
        *m_lightEventMetaDataDescriptors[mdInd++]);
    crossingPointVal->setValue(ts);
    mdArray.add(crossingPointVal);

    MetaDataValue* directionVal = new MetaDataValue(
        *m_lightEventMetaDataDescriptors[mdInd++]);
    directionVal->setValue(static_cast<juce::uint8>(val));
    mdArray.add(directionVal);

    // MetaDataValue* crossingPointVal = new MetaDataValue(
    //     *m_thetaEventMetaDataDescriptors[mdInd++]);
    // crossingPointVal->setValue(ts);
    // mdArray.add(crossingPointVal);

    // MetaDataValue* threshVal = new MetaDataValue(
    //     *m_thetaEventMetaDataDescriptors[mdInd++]);
    // threshVal->setValue(m_threshold);
    // mdArray.add(threshVal);

    // MetaDataValue* directionVal = new MetaDataValue(
    //     *m_thetaEventMetaDataDescriptors[mdInd++]);
    // directionVal->setValue(val);
    // mdArray.add(directionVal);

    uint8 ttlData;
    ttlData = val << m_lightEventChan;

    TTLEventPtr event = TTLEvent::createTTLEvent(
        m_lightEventChannelPtr,
        // m_thetaEventChannelPtr,
        ts, &ttlData, sizeof(uint8), mdArray, m_lightEventChan
        );
    addEvent(m_lightEventChannelPtr, event, tsInBuffer);
    // addEvent(m_thetaEventChannelPtr, event, tsInBuffer);

    if (val) {
        // calculate time to off
        // int sampRate = upstreamDataChannelArray[m_inputChan]->getSampleRate();
        int sampRate = getDataChannel(0)->getSampleRate();
        float dur = m_lightDur * sampRate;

        m_tsLightOn = ts;
        m_tsLightOff = ts + static_cast<juce::int64>(dur);
        m_isLightOn = true;
        m_tsDelayLight = -1;

    } else {
        m_isLightOn = false;
    }
}

// all new values should be validated before this function is called!
void ThetaDetector::setParameter(int parameterIndex, float newValue) {
    switch (parameterIndex) {
    case pInputChan:
        if (getNumInputs() > newValue)
            m_inputChan = static_cast<int>(newValue);
        break;

	case pThetaLow:
		m_thetaLow = newValue;
        m_theta_estimator.setThetaLow(m_thetaLow);
		setFilterParams();
		break;

	case pThetaHigh:
		m_thetaHigh = newValue;
        m_theta_estimator.setThetaHigh(m_thetaHigh);
		setFilterParams();
		break;

	case pRefLow:
		m_refLow = newValue;
		setFilterParams();
		break;

	case pRefHigh:
		m_refHigh = newValue;
		setFilterParams();
		break;

	case pThreshold:
		m_threshold = newValue;
        m_thetaMonitor.setThreshold(m_threshold);
		break;

	case pDelayEnabled:
		m_isDelayEnabled = newValue;
		break;

    case pLightGateChan:
        m_lightGateChan = newValue;
        break;
    }
}

bool ThetaDetector::disable() {
    return true;
}
