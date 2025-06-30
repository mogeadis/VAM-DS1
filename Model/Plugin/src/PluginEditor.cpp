#include "PluginEditor.h"

// ====================================================================================================
// Class: PluginEditor
// ====================================================================================================

// [Public] Constructor & Destructor Methods
PluginEditor::PluginEditor(PluginProcessor& processor) : AudioProcessorEditor(&processor),processor(processor)
{
    juce::ignoreUnused(processor);
    setSize(400,300);

    // Distortion Slider
    addAndMakeVisible(this->distortion_slider);
    this->distortion_slider.onValueChange = [this]
    {
        this->processor.distortion = static_cast<float>(this->distortion_slider.getValue());
    };
    this->distortion_slider.setRange(0,1,0.01);
    this->distortion_slider.setValue(0.5);
    this->distortion_slider.setSliderStyle(juce::Slider::LinearVertical);
    this->distortion_slider.setTextBoxStyle(juce::Slider::NoTextBox,false,0,0);
    this->distortion_slider.setPopupDisplayEnabled(true,false,this);

    // Tone Slider
    addAndMakeVisible(this->tone_slider);
    this->tone_slider.onValueChange = [this]
    {
        this->processor.tone = static_cast<float>(this->tone_slider.getValue());
    };
    this->tone_slider.setRange(0,1,0.01);
    this->tone_slider.setValue(0.5);
    this->tone_slider.setSliderStyle(juce::Slider::LinearVertical);
    this->tone_slider.setTextBoxStyle(juce::Slider::NoTextBox,false,0,0);
    this->tone_slider.setPopupDisplayEnabled(true,false,this);

    // Level Slider
    addAndMakeVisible(this->level_slider);
    this->level_slider.onValueChange = [this]
    {
        this->processor.level = static_cast<float>(this->level_slider.getValue());
    };
    this->level_slider.setRange(0,1,0.01);
    this->level_slider.setValue(0.5);
    this->level_slider.setSliderStyle(juce::Slider::LinearVertical);
    this->level_slider.setTextBoxStyle(juce::Slider::NoTextBox,false,0,0);
    this->level_slider.setPopupDisplayEnabled(true,false,this);

    // Bypass Button
    addAndMakeVisible(this->bypass_button);
    this->bypass_button.setButtonText("Bypass");
    this->bypass_button.setClickingTogglesState(true);
    this->bypass_button.onClick = [this]
    {
        if(this->bypass_button.getToggleState())
        {
            this->processor.bypass = true;
        }
        else
        {
            this->processor.bypass = false;
        }
    };
}
PluginEditor::~PluginEditor(){}

// [Public] JUCE Default Methods
void PluginEditor::paint(juce::Graphics& graphics)
{
    graphics.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));
    graphics.setColour(juce::Colours::white);
    graphics.setFont(15.0f);
}
void PluginEditor::resized()
{
    int x,y,width,height;
    this->window_width = static_cast<float>(getWidth());
    this->window_height = static_cast<float>(getHeight());

    // Distortion Slider
    width = 30;
    height = 200;
    x = int(this->window_width/4);
    y = int(this->window_height/2);
    this->distortion_slider.setSize(width,height);
    this->distortion_slider.setCentrePosition(x,y);

    // Tone Slider
    width = 30;
    height = 200;
    x = int(this->window_width/2);
    y = int(this->window_height/2);
    this->tone_slider.setSize(width,height);
    this->tone_slider.setCentrePosition(x,y);

    // Level Slider
    width = 30;
    height = 200;
    x = int(this->window_width*3/4);
    y = int(this->window_height/2);
    this->level_slider.setSize(width,height);
    this->level_slider.setCentrePosition(x,y);

    // Bypass Button
    width = 100;
    height = 30;
    x = int(this->window_width/2);
    y = int(this->window_height*5/6);
    this->bypass_button.setSize(width,height);
    this->bypass_button.setCentrePosition(x,y);
}