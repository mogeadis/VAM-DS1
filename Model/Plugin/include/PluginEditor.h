#pragma once

#include <JuceHeader.h>

#include "PluginProcessor.h"

// ====================================================================================================
// Class: PluginEditor
// ====================================================================================================
class PluginEditor final : public juce::AudioProcessorEditor
{
    public:
    // Constructor & Destructor Methods
    PluginEditor(PluginProcessor& processor);
    ~PluginEditor() override;

    // JUCE Default Methods
    void paint(juce::Graphics& graphics) override;
    void resized() override;

    private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PluginEditor)

    // Member Variables
    float window_width;
    float window_height;
    PluginProcessor& processor;
    juce::Slider distortion_slider;
    juce::Slider tone_slider;
    juce::Slider level_slider;
    juce::TextButton bypass_button;
};