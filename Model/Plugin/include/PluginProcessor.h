#pragma once
#pragma warning(disable : 5030)

#include <JuceHeader.h>

#include <atomic>

#include "bossds1.h"

// ====================================================================================================
// Class: PluginProcessor
// ====================================================================================================
class PluginProcessor final : public juce::AudioProcessor
{
    public:
    using AudioProcessor::processBlock;

    // Member Variables
    std::unique_ptr<BossDS1> effect;
    std::atomic<float> distortion;
    std::atomic<float> tone;
    std::atomic<float> level;
    std::atomic<bool> bypass;

    // Constructor & Destructor Methods
    PluginProcessor();
    ~PluginProcessor() override;

    // JUCE Default Methods
    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    void prepareToPlay(double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
    const juce::String getName() const override;
    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram(int index) override;
    const juce::String getProgramName(int index) override;
    void changeProgramName(int index, const juce::String& newName) override;
    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;
    bool hasEditor() const override;
    juce::AudioProcessorEditor* createEditor() override;
    void getStateInformation(juce::MemoryBlock& destData) override;
    void setStateInformation(const void* data, int sizeInBytes) override;

    private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PluginProcessor)

    // Member Variables
    int buffer_samples;
    int input_channels;
    int output_channels;
};
