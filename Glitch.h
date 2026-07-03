#pragma once

#include "vessl/vessl.h"

using count_t = uint32_t;

struct FreezeSettings
{
  // used to determine how long the frozen section of audio should be.
  float clockRatio;
  // how many clock ticks should occur before resetting the read LFO when not frozen,
  // in order to keep it in sync with the clock.
  count_t readResetCount;
  // param value at which to choose this setting
  float paramThresh;
};

static const FreezeSettings FREEZE_SETTINGS[] = {
  { 0.25f, 1, 0.0f  },
  { 0.5f,  1, 0.02f },
  { 1.0f,  1, 0.06f },
  { 1.5f,  3, 0.20f },
  { 2.0f,  2, 0.4f  },
  { 3.0f,  3, 0.6f  },
  { 4.0f,  4, 0.7f  },
  { 6.0f,  6, 0.85f },
  { 8.0f,  8, 0.95f },
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
             , protected vessl::plist<7>
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
    vessl::analog_p play_rate;
    vessl::binary_p freeze;
    vessl::binary_p glitchEnabled;
  } params_;
  BufferType freeze_buffer_;
  Freeze freeze_proc_;
    
  float sample_rate_;
  float glitch_lfo_;
  float glitch_rand_;
  float freeze_rate_;
  count_t freeze_settings_idx_;
  count_t glitch_settings_idx_;
  count_t freeze_counter_;
  count_t glitch_counter_;
  count_t samples_since_last_tap_;
  
  BitCrush crush_proc_;
  
  BufferType process_buffer_;

  Array follower_window_;
  EnvelopeFollower envelope_follower_;
  Array input_envelope_;
  
public:
  Glitch(float sample_rate, vessl::size_t block_size)
  : clockable(sample_rate, static_cast<uint32_t>(block_size), FREEZE_BUFFER_SIZE)
  , freeze_buffer_(new GlitchSampleType[FREEZE_BUFFER_SIZE], FREEZE_BUFFER_SIZE)
  , freeze_proc_(freeze_buffer_, sample_rate)
  , sample_rate_(sample_rate)
  , glitch_lfo_(0)
  , glitch_rand_(0)
  , freeze_rate_(0)
  , freeze_settings_idx_(0)
  , glitch_settings_idx_(0)
  , freeze_counter_(0), glitch_counter_(0)
  , samples_since_last_tap_(FREEZE_BUFFER_SIZE)
  , crush_proc_(sample_rate, sample_rate)
  , process_buffer_(new GlitchSampleType[block_size], block_size)
  , follower_window_(new float[block_size * 8], block_size * 8) // NOLINT(bugprone-implicit-widening-of-multiplication-result)
  , envelope_follower_(follower_window_, sample_rate, 0.001f)
  , input_envelope_(new float[block_size], block_size)
  {
  }

  ~Glitch() override
  {
    delete[] input_envelope_.data();
    delete[] follower_window_.data();
    delete[] freeze_buffer_.data();
    delete[] process_buffer_.data();
  }

  using clockable::clock;

  [[nodiscard]] parameter repeats() const { return params_.repeats("repeats", 'r');  }
  [[nodiscard]] parameter crush() const { return params_.crush("crush", 'c'); }
  [[nodiscard]] parameter glitch() const { return params_.glitch("glitch", 'g'); }
  [[nodiscard]] parameter glitching() const { return params_.glitchEnabled("glich enabled", 'e'); }
  [[nodiscard]] parameter shape() const { return params_.shape("shape", 's'); }
  [[nodiscard]] parameter freeze() const { return params_.freeze("freeze", 'f'); }
  [[nodiscard]] parameter play_rate() const { return params_.play_rate("play rate", 'p'); }
  [[nodiscard]] float freeze_phase() const { return freeze_proc_.phase(); }
  [[nodiscard]] float envelope() const { return input_envelope_[0]; }
  [[nodiscard]] float glitch_rand() const { return glitch_rand_; }

  void process(vessl::array<GlitchSampleType> input, vessl::array<GlitchSampleType> output) override
  {
    vessl::size_t size = input.size();
    clockable::tick(size);
    
    float smooth_freeze = repeats();
    for (freeze_settings_idx_ = 0; freeze_settings_idx_ < FREEZE_SETTINGS_COUNT - 1; freeze_settings_idx_++)
    {
      if (smooth_freeze >= FREEZE_SETTINGS[freeze_settings_idx_].paramThresh
        && smooth_freeze < FREEZE_SETTINGS[freeze_settings_idx_+1].paramThresh)
      {
        break;
      }
    }
    
    float new_freeze_length = freeze_size(freeze_settings_idx_);
    //float newReadSpeed = freeze_speed(freeze_settings_idx_);
    float new_read_speed = play_rate().read_analog();
    
    // smooth size and speed changes when not clocked
    bool clocked = samples_since_last_tap_ < FREEZE_BUFFER_SIZE;
    if (!clocked)
    {
      if (freeze_settings_idx_ < FREEZE_SETTINGS_COUNT - 1)
      {
        float p0 = FREEZE_SETTINGS[freeze_settings_idx_].paramThresh;
        float p1 = FREEZE_SETTINGS[freeze_settings_idx_+1].paramThresh;
        float t = (smooth_freeze - p0) / (p1 - p0);
        float d1 = freeze_size(freeze_settings_idx_ + 1);
        new_freeze_length = new_freeze_length + (d1 - new_freeze_length)*t;
      }
    }
    
    freeze_proc_.duration() = new_freeze_length;
    freeze_proc_.rate() = new_read_speed;
    freeze_proc_.enabled() = freeze().read_binary();
    
    float sr = sample_rate_;
    float crush_param = crush();
    float bits = crush_param > 0.001f ? (16.f - crush_param*12.0f) : 24;
    float rate = crush_param > 0.001f ? sr * 0.25f + crush_param*(100 - sr * 0.25f) : sr;
    crush_proc_.depth() = bits;
    crush_proc_.rate() = rate;
    
    auto input_reader = input.make_reader();
    auto procw = process_buffer_.make_writer();
    auto iew = input_envelope_.make_writer();
    while(input_reader)
    {
      GlitchSampleType sample = input_reader.read();
      procw << sample;
      iew << sample.to_mono().value();
    }
    envelope_follower_.process(input_envelope_, input_envelope_);
    
    //can't use output as a process buffer because we need the dry input again for the shape stage.
    if (clocked)  // NOLINT(bugprone-branch-clone)
    {
      freeze_proc_.process<vessl::time::mode::fade>(process_buffer_, process_buffer_);
    }
    else
    {
      freeze_proc_.process<vessl::time::mode::slew>(process_buffer_, process_buffer_);
    }
    
    crush_proc_.process(process_buffer_, process_buffer_);
    
    float glitch_param = glitch();
    glitch_settings_idx_ = static_cast<int>((1.f - glitch_param) * GLITCH_SETTINGS_COUNT);
    float glitch_speed = 1.0f / (glitch_size(glitch_settings_idx_) * GLITCH_LFO_DIV);
    float glitch_prob = glitch_param < 0.001f ? 0 : 0.1f + 0.4f*glitch_param;
    if (glitch_prob == 0)
    {
      params_.glitchEnabled.value = false;
    }
    for (count_t i = 0; i < size; ++i)
    {
      if (step_glitch_lfo(glitch_speed))
      {
        glitch_rand_ = vessl::math::random::range<float>(0.f, 1.f);
        if (glitch_rand_ < glitch_prob)
        {
          params_.glitchEnabled.value = !params_.glitchEnabled.value;
        }
        //params.glitchEnabled.value = glitchRand < glitch_prob;
      }
    
      if (params_.glitchEnabled.value)
      {
        vessl::size_t d = i+1;
        GlitchSampleType f = freeze_proc_.buffer().read(d);
        GlitchSampleType& pf = process_buffer_[i];
        pf.left() = glitch(pf.left(), f.left());
        pf.right() = glitch(pf.right(), f.right());
      }
    }
    
    float shape_param = shape();
    float shape_wet = shape_param;
    float shape_dry = 1.0f - shape_wet;
    float size_f = static_cast<float>(size);
    input_reader.reset();
    for (count_t i = 0; i < size; ++i)
    {
      const float shape_scale = input_envelope_[i]*size_f*(10.0f + 90.0f*shape_param);
      const float dry_idx = static_cast<float>(i);
      // treat the process buffer like a wave table and use the dry input as phase, modulated by the envelope follower,
      // using shapeParam both for dry/wet mix and scaling of the envelope value.
      GlitchSampleType in = input_reader.read();
      const float read_l = shape_dry*dry_idx + shape_wet*vessl::math::constrain(shape_scale*in.left(), -size_f, size_f);
      const float read_r = shape_dry*dry_idx + shape_wet*vessl::math::constrain(shape_scale*in.right(), -size_f, size_f);
      output[i] = {interpolatedReadAt(process_buffer_, read_l).left(), interpolatedReadAt(process_buffer_, read_r).right()};
    }
    
    if (samples_since_last_tap_ < FREEZE_BUFFER_SIZE)
    {
      samples_since_last_tap_ += size;
    }
  }

protected:
  [[nodiscard]] parameter element_at(vessl::size_t index) const override
  {
    parameter p[num] = { repeats(), crush(), glitch(), glitching(), shape(), freeze(), play_rate() };
    return p[index];
  }
  
  void tock(period_t sample_delay) override
  {
    samples_since_last_tap_ = 0;
      
    // reset readLfo based on the counter for our current setting
    if (++freeze_counter_ >= FREEZE_SETTINGS[freeze_settings_idx_].readResetCount)
    {
      freeze_proc_.reset();
      freeze_counter_ = 0;
    }

    // // we use one instead of zero because our logic in process
    // // is checking for the flip from 1 to 0 to generate a new random value.
    if (++glitch_counter_ >= GLITCH_SETTINGS[glitch_settings_idx_].lfoResetCount*GLITCH_LFO_DIV)
    {
      glitch_lfo_ = 1;
      glitch_counter_ = 0;
    }
  }

private:
  VESSL_INLINE bool step_glitch_lfo(const float speed)
  {
    glitch_lfo_ = glitch_lfo_ + speed;
    if (glitch_lfo_ >= 1)
    {
      glitch_lfo_ -= 1;
      return true;
    }
    if (glitch_lfo_ < 0)
    {
      glitch_lfo_ += 1;
      return true;
    }
    return false;
  }

  [[nodiscard]] VESSL_INLINE float freeze_size(const count_t idx) const
  {
    return period() * FREEZE_SETTINGS[idx].clockRatio;
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