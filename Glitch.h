#pragma once

#include "vessl/vessl.h"

using count_t = uint32_t;

struct FreezeSettings
{
  // used to determine how long the frozen section of audio should be.
  float clockRatio;
  // these are the speeds at which the frozen audio should be played back.
  float playbackSpeed;
  // how many clock ticks should occur before resetting the read LFO when not frozen,
  // in order to keep it in sync with the clock.
  count_t readResetCount;
  // param value at which to choose this setting
  float paramThresh;
};

static const FreezeSettings FREEZE_SETTINGS[] = {
  { 2.0f,     4.0f, 1, 0.0f  },
  { 2.0f,     3.0f, 2, 0.02f },
  { 2.0f,     2.0f, 1, 0.06f },
  { 4.0f/3.0f,1.0f, 4, 0.20f },
  { 2.0f,     1.0f, 2, 0.4f  },
  { 3.0f,     1.0f, 3, 0.6f  },
  { 4.0f,     1.0f, 4, 0.7f  },
  { 6.0f,     1.0f, 6, 0.85f },
  { 8.0f,     1.0f, 8, 0.95f },
};
static constexpr count_t FREEZE_SETTINGS_COUNT = sizeof(FREEZE_SETTINGS) / sizeof(FreezeSettings);

struct GlitchSettings
{
  float clockRatio;
  count_t lfoResetCount;
};

static constexpr GlitchSettings GLITCH_SETTINGS[] = {
  { 1.0f / 32, 1 },
  { 1.0f / 24, 1 },
  { 1.0f / 16, 1 },
  { 1.0f / 12, 1 },
  { 1.0f / 8, 1 },
  { 1.0f / 6, 1 },
  { 1.0f / 4, 1 },
  { 1.0f / 3, 1 },
  { 1.0f / 2, 1 },
  { 1, 1 },
};
static constexpr count_t GLITCH_SETTINGS_COUNT = sizeof(GLITCH_SETTINGS) / sizeof(GlitchSettings);

static constexpr count_t GLITCH_LFO_DIV = 4;

using GlitchSampleType = vessl::sample::type<float>::stereo;
using BufferType = vessl::array<GlitchSampleType>;
using BitCrush = vessl::processors::bitcrush<GlitchSampleType, 24>;
using Freeze = vessl::processors::freeze<GlitchSampleType>;
using EnvelopeFollower = vessl::processors::follow<float>;
using Array = vessl::array<float>;

template<uint32_t FREEZE_BUFFER_SIZE>
class Glitch : public vessl::unit_processor<GlitchSampleType>
             , public vessl::time::clockable
             , protected vessl::plist<6>
{
public:
  using parameter = vessl::parameter;
  const parameter_list& parameters() const override { return *this; }
 
private:
  struct
  {
    vessl::analog_p repeats;
    vessl::analog_p crush;
    vessl::analog_p glitch;
    vessl::analog_p shape;
    vessl::binary_p freeze;
    vessl::binary_p glitchEnabled;
  } params;
  BufferType freezeBuffer;
  Freeze freezeProc;
    
  float sampleRate;
  float glitchLfo;
  float glitchRand;
  count_t freezeSettingsIdx;
  count_t glitchSettingsIdx;
  count_t freezeCounter;
  count_t glitchCounter;
  count_t samplesSinceLastTap;
  
  BitCrush crushProc;
  
  BufferType processBuffer;

  Array followerWindow;
  EnvelopeFollower envelopeFollower;
  Array inputEnvelope;
  
public:
  Glitch(float sampleRate, vessl::size_t blockSize) 
  : clockable(sampleRate, static_cast<uint32_t>(blockSize), FREEZE_BUFFER_SIZE)
  , freezeBuffer(new GlitchSampleType[FREEZE_BUFFER_SIZE], FREEZE_BUFFER_SIZE)
  , freezeProc(freezeBuffer, sampleRate)
  , sampleRate(sampleRate)
  , glitchLfo(0), glitchRand(0), freezeSettingsIdx(0), glitchSettingsIdx(0)
  , freezeCounter(0), glitchCounter(0)
  , samplesSinceLastTap(FREEZE_BUFFER_SIZE)
  , crushProc(sampleRate, sampleRate)
  , processBuffer(new GlitchSampleType[blockSize], blockSize)
  , followerWindow(new float[blockSize*8], blockSize*8)  // NOLINT(bugprone-implicit-widening-of-multiplication-result)
  , envelopeFollower(followerWindow, sampleRate, 0.001f)
  , inputEnvelope(new float[blockSize], blockSize)
  {
  }
  
  ~Glitch() override
  {
    delete[] inputEnvelope.data();
    delete[] followerWindow.data();
    delete[] freezeBuffer.data();
    delete[] processBuffer.data();
  }

  using clockable::clock;

  [[nodiscard]] parameter repeats() const { return params.repeats("repeats", 'r');  }
  [[nodiscard]] parameter crush() const { return params.crush("crush", 'c'); }
  [[nodiscard]] parameter glitch() const { return params.glitch("glitch", 'g'); }
  [[nodiscard]] parameter glitching() const { return params.glitchEnabled("glich enabled", 'e'); }
  [[nodiscard]] parameter shape() const { return params.shape("shape", 's'); }
  [[nodiscard]] parameter freeze() const { return params.freeze("freeze", 'f'); }
  [[nodiscard]] float freeze_phase() const { return freezeProc.phase(); }
  [[nodiscard]] float envelope() const { return inputEnvelope[0]; }
  [[nodiscard]] float glitch_rand() const { return glitchRand; }

  void process(vessl::array<GlitchSampleType> input, vessl::array<GlitchSampleType> output) override
  {
    vessl::size_t size = input.size();
    clockable::tick(size);
    
    float smoothFreeze = repeats();
    for (freezeSettingsIdx = 0; freezeSettingsIdx < FREEZE_SETTINGS_COUNT - 1; freezeSettingsIdx++)
    {
      if (smoothFreeze >= FREEZE_SETTINGS[freezeSettingsIdx].paramThresh
        && smoothFreeze < FREEZE_SETTINGS[freezeSettingsIdx+1].paramThresh)
      {
        break;
      }
    }
    
    float newFreezeLength = freeze_size(freezeSettingsIdx);
    float newReadSpeed = freeze_speed(freezeSettingsIdx);
    
    // smooth size and speed changes when not clocked
    bool clocked = samplesSinceLastTap < FREEZE_BUFFER_SIZE;
    if (!clocked)
    {
      if (freezeSettingsIdx < FREEZE_SETTINGS_COUNT - 1)
      {
        float p0 = FREEZE_SETTINGS[freezeSettingsIdx].paramThresh;
        float p1 = FREEZE_SETTINGS[freezeSettingsIdx+1].paramThresh;
        float t = (smoothFreeze - p0) / (p1 - p0);
        float d1 = freeze_size(freezeSettingsIdx + 1);
        newFreezeLength = newFreezeLength + (d1 - newFreezeLength)*t;
        newReadSpeed = newReadSpeed + (freeze_speed(freezeSettingsIdx + 1) - newReadSpeed)*t;
      }
    }
    
    freezeProc.duration() = newFreezeLength;
    freezeProc.rate() = newReadSpeed;
    freezeProc.enabled() = freeze().read_binary();
    
    float sr = sampleRate;
    float crushParam = crush();
    float bits = crushParam > 0.001f ? (16.f - crushParam*12.0f) : 24;
    float rate = crushParam > 0.001f ? sr * 0.25f + crushParam*(100 - sr * 0.25f) : sr;
    crushProc.depth() = bits;
    crushProc.rate() = rate;
    
    auto inputReader = input.make_reader();
    auto procw = processBuffer.make_writer();
    auto iew = inputEnvelope.make_writer();
    while(inputReader)
    {
      GlitchSampleType sample = inputReader.read();
      procw << sample;
      iew << sample.to_mono().value();
    }
    envelopeFollower.process(inputEnvelope, inputEnvelope);
    
    //can't use output as a process buffer because we need the dry input again for the shape stage.
    if (clocked)  // NOLINT(bugprone-branch-clone)
    {
      freezeProc.process<vessl::time::mode::fade>(processBuffer, processBuffer);
    }
    else
    {
      freezeProc.process<vessl::time::mode::slew>(processBuffer, processBuffer);
    }
    
    crushProc.process(processBuffer, processBuffer);
    
    float glitch_param = glitch();
    glitchSettingsIdx = static_cast<int>((1.f - glitch_param) * GLITCH_SETTINGS_COUNT);
    float glitch_speed = 1.0f / (glitch_size(glitchSettingsIdx) * GLITCH_LFO_DIV);
    float glitch_prob = glitch_param < 0.001f ? 0 : 0.1f + 0.4f*glitch_param;
    if (glitch_prob == 0)
    {
      params.glitchEnabled.value = false;
    }
    for (count_t i = 0; i < size; ++i)
    {
      if (step_glitch_lfo(glitch_speed))
      {
        glitchRand = vessl::math::random::range<float>(0.f, 1.f);
        if (glitchRand < glitch_prob)
        {
          params.glitchEnabled.value = !params.glitchEnabled.value;
        }
        //params.glitchEnabled.value = glitchRand < glitch_prob;
      }
    
      if (params.glitchEnabled.value)
      {
        vessl::size_t d = i+1;
        GlitchSampleType f = freezeProc.buffer().read(d);
        GlitchSampleType& pf = processBuffer[i];
        pf.left() = glitch(pf.left(), f.left());
        pf.right() = glitch(pf.right(), f.right());
      }
    }
    
    float shapeParam = shape();
    float shapeWet = shapeParam;
    float shapeDry = 1.0f - shapeWet;
    float fSize = static_cast<float>(size);
    inputReader.reset();
    for (count_t i = 0; i < size; ++i)
    {
      const float shapeScale = inputEnvelope[i]*fSize*(10.0f + 90.0f*shapeParam);
      const float dryIdx = static_cast<float>(i);
      // treat the process buffer like a wave table and use the dry input as phase, modulated by the envelope follower,
      // using shapeParam both for dry/wet mix and scaling of the envelope value.
      GlitchSampleType in = inputReader.read();
      const float readL = shapeDry*dryIdx + shapeWet*vessl::math::constrain(shapeScale*in.left(), -fSize, fSize);
      const float readR = shapeDry*dryIdx + shapeWet*vessl::math::constrain(shapeScale*in.right(), -fSize, fSize);
      output[i] = GlitchSampleType(interpolatedReadAt(processBuffer, readL).left(), interpolatedReadAt(processBuffer, readR).right());
    }
    
    if (samplesSinceLastTap < FREEZE_BUFFER_SIZE)
    {
      samplesSinceLastTap += size;
    }
  }

protected:
  [[nodiscard]] parameter element_at(vessl::size_t index) const override
  {
    parameter p[num] = { repeats(), crush(), glitch(), glitching(), shape(), freeze() };
    return p[index];
  }
  
  void tock(period_t sample_delay) override
  {
    samplesSinceLastTap = 0;
      
    // reset readLfo based on the counter for our current setting
    if (++freezeCounter >= FREEZE_SETTINGS[freezeSettingsIdx].readResetCount)
    {
      freezeProc.reset();
      freezeCounter = 0;
    }

    // // we use one instead of zero because our logic in process
    // // is checking for the flip from 1 to 0 to generate a new random value.
    if (++glitchCounter >= GLITCH_SETTINGS[glitchSettingsIdx].lfoResetCount*GLITCH_LFO_DIV)
    {
      glitchLfo = 1;
      glitchCounter = 0;
    }
  }

private:
  VESSL_INLINE bool step_glitch_lfo(const float speed)
  {
    glitchLfo = glitchLfo + speed;
    if (glitchLfo >= 1)
    {
      glitchLfo -= 1;
      return true;
    }
    if (glitchLfo < 0)
    {
      glitchLfo += 1;
      return true;
    }
    return false;
  }

  [[nodiscard]] VESSL_INLINE float freeze_size(const count_t idx) const
  {
    return period() * FREEZE_SETTINGS[idx].clockRatio;
  }

  VESSL_INLINE static float freeze_speed(const count_t idx)
  {
    return FREEZE_SETTINGS[idx].playbackSpeed;
  }

  [[nodiscard]] VESSL_INLINE float glitch_size(const count_t idx) const
  {
    return period() * GLITCH_SETTINGS[idx].clockRatio;
  }

  VESSL_INLINE static float glitch(const float a, const float b)
  {
    vessl::q31_t glitched = vessl::cast<vessl::q31_t>(a) ^ vessl::cast<vessl::q31_t>(b);
    return vessl::cast<float>(glitched);
  }
  
  VESSL_INLINE static GlitchSampleType interpolatedReadAt(vessl::array<GlitchSampleType> buffer, float index)
  {
    // index can be negative, we ensure it is positive.
    index += static_cast<float>(buffer.size());
    count_t idx = static_cast<count_t>(index);
    const GlitchSampleType& low = buffer[idx%buffer.size()];
    const GlitchSampleType& high = buffer[(idx + 1)%buffer.size()];
    float frac = index - static_cast<float>(idx);
    return low + frac * (high - low);
  }
  
  GlitchSampleType process(const GlitchSampleType& in) override { return in; }
};