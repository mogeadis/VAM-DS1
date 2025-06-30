#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PluginProcessor::PluginProcessor(): AudioProcessor
(BusesProperties()
#if ! JucePlugin_IsMidiEffect
    #if ! JucePlugin_IsSynth
        .withInput("Input",juce::AudioChannelSet::stereo(),true)
    #endif
    .withOutput("Output",juce::AudioChannelSet::stereo(),true)
#endif
)
{
    this->distortion = 0.5f;
    this->tone = 0.5f;
    this->level = 0.5f;
    this->bypass = false;
    this->input_channels  = getTotalNumInputChannels();
    this->output_channels = getTotalNumOutputChannels();
}
PluginProcessor::~PluginProcessor(){}

// [Public] JUCE Default Methods
void PluginProcessor::processBlock(juce::AudioBuffer<float>& buffer,juce::MidiBuffer& midiMessages)
{
    juce::ignoreUnused(midiMessages);
    juce::ScopedNoDenormals noDenormals;
    this->effect->setDistortion(this->distortion);
    this->effect->setTone(this->tone);
    this->effect->setLevel(this->level);
    for(int n = 0; n < this->buffer_samples; n++)
    {
        float Vin = 0;
        for (int channel = 0; channel < input_channels; channel++)
        {
            Vin += buffer.getReadPointer(channel)[n];
        }
        Vin = Vin/input_channels;
        float Vout = this->effect->processSample(Vin);
        for (int channel = 0; channel < output_channels; channel++)
        {
            if(this->bypass)
            {
                buffer.getWritePointer(channel)[n] = Vin;
            }
            else
            {
                buffer.getWritePointer(channel)[n] = Vout;
            }
        }
    }
}
void PluginProcessor::prepareToPlay(double sampleRate,int samplesPerBlock)
{
    this->buffer_samples = samplesPerBlock;
    this->effect.reset(new BossDS1(static_cast<float>(sampleRate),this->distortion,this->tone,this->level));
}
void PluginProcessor::releaseResources(){}
const juce::String PluginProcessor::getName() const
{
    return JucePlugin_Name;
}
bool PluginProcessor::acceptsMidi() const
{
    #if JucePlugin_WantsMidiInput
        return true;
    #else
        return false;
    #endif
}
bool PluginProcessor::producesMidi() const
{
    #if JucePlugin_ProducesMidiOutput
        return true;
    #else
        return false;
    #endif
}
bool PluginProcessor::isMidiEffect() const
{
    #if JucePlugin_IsMidiEffect
        return true;
    #else
        return false;
    #endif
}
double PluginProcessor::getTailLengthSeconds() const
{
    return 0.0;
}
int PluginProcessor::getNumPrograms()
{
    return 1;
}
int PluginProcessor::getCurrentProgram()
{
    return 0;
}
void PluginProcessor::setCurrentProgram (int index)
{
    juce::ignoreUnused(index);
}
const juce::String PluginProcessor::getProgramName(int index)
{
    juce::ignoreUnused(index);
    return {};
}
void PluginProcessor::changeProgramName(int index,const juce::String& newName)
{
    juce::ignoreUnused(index,newName);
}
bool PluginProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    #if JucePlugin_IsMidiEffect
        juce::ignoreUnused(layouts);
        return true;
    #else
        if(layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono() && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
            return false;
        #if ! JucePlugin_IsSynth
            if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
                return false;
        #endif
        return true;
    #endif
}
bool PluginProcessor::hasEditor() const
{
    return true;
}
juce::AudioProcessorEditor* PluginProcessor::createEditor()
{
    return new PluginEditor(*this);
}
void PluginProcessor::getStateInformation(juce::MemoryBlock& destData)
{
    juce::ignoreUnused(destData);
}
void PluginProcessor::setStateInformation(const void* data, int sizeInBytes)
{
    juce::ignoreUnused(data,sizeInBytes);
}
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new PluginProcessor();
}