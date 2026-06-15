#pragma once

#include "vessl/vessl.h"

template<typename T, unsigned ChannelCount, unsigned MaxGrains>
class Granulator : public vessl::unit_processor<vessl::sample::frame<T, ChannelCount>>
                 , public vessl::generator<vessl::sample::frame<T, ChannelCount>>
                 , protected vessl::plist<5>
{
public:
  using Parameter = vessl::parameter;
  using SampleType = vessl::sample::frame<T, ChannelCount>;
  using RecordSampleType = typename vessl::sample::type<T>::mono;
  
  [[nodiscard]] const parameter_list& parameters() const override { return *this; }
  
  // length of a grain in seconds (duration_t)
  [[nodiscard]] Parameter grain_duration() const { return params_.grain_duration("g.dur", 'd'); }
  // playback speed of grain, clamped to [-2,2]
  [[nodiscard]] Parameter grain_speed() const { return params_.grain_speed("g.spd", 's'); }
  // how far back in time the grain start should be moved (duration_t)
  [[nodiscard]] Parameter grain_offset() const { return params_.grain_offset("g.off", 'o'); }
  // how often grains are started (duration_t)
  [[nodiscard]] Parameter grain_rate() const { return params_.grain_rate("g.rate", 'r'); }
  // how a grain is panned in the multi-channel field [-1,1]
  [[nodiscard]] Parameter grain_pan() const { return params_.grain_pan("g.pan", 'p'); }
  
  // starts a new grain
  VESSL_INLINE void trigger(float sample_delay = 0)
  {
    if (active_grain_count_ < MaxGrains)
    {
      Grain& grain = grains_[active_grain_count_++];
      grain.speed = vessl::math::constrain(params_.grain_speed.value, 0.5f, 2.f);
      grain.size = params_.grain_duration.value.samples * grain.speed;
      grain.start = record_buffer_.get_write_index()
                  - grain.size 
                  - params_.grain_offset.value.samples
                  // make sure we're working with positive indices
                  + record_buffer_.size();
      grain.pan = vessl::math::constrain(params_.grain_pan.value, -1.f, 1.f);
      
      // for now
      float env = 0.5f;
      float next_attack = vessl::math::constrain(env, 0.01f, 0.99f);
      float next_decay = 1.0f - next_attack;
      grain.decay_start = next_attack * grain.size;
      grain.attack_mult = 1.0f / (next_attack * grain.size);
      grain.decay_mult = 1.0f / (next_decay * grain.size);
      
      grain.ramp = -sample_delay;
      
      grain_triggered_ = true;
    }
  }
    
  // should only be called immediately after process/generate
  [[nodiscard]] bool started_grain() const { return grain_triggered_; }
  
  [[nodiscard]] int active_grain_count() const { return active_grain_count_;}

  [[nodiscard]] VESSL_INLINE SampleType process(const SampleType &in) override
  {
    grain_triggered_ = false;
    record_buffer_.write(in.to_mono());
    return generate();
  }
  
  VESSL_INLINE void process(const vessl::array<SampleType>& in, vessl::array<SampleType> out)
  {
    grain_triggered_ = false;
    auto rin = in.make_reader();
    //auto wout = out.make_writer();
    while (rin)
    {
      SampleType s = rin.read();
      record_buffer_.write(s.to_mono());
      //wout << generate();
    }
    
    generate(out);
  }
  
  [[nodiscard]] VESSL_INLINE SampleType generate() override
  {
    grain_rate_phasor_++;
    
    if (grain_rate_phasor_ >= params_.grain_rate.value.samples)
    {
      trigger();
      grain_rate_phasor_ = 0;
    }
    
    SampleType accum = SampleType(0);
    SampleType samp;
    RecordSampleType* buffer = record_buffer_.data();
    for (int i = active_grain_count_ - 1; i >= 0; i--)
    {
      Grain& grain = grains_[i];
      if (grain.ramp >= 0)
      {
        float pos = grain.start + grain.ramp;
        float env = grain.ramp < grain.decay_start 
                  ? grain.ramp * grain.attack_mult 
                  : (grain.size - grain.ramp) * grain.decay_mult;
        int i = static_cast<int>(pos);
        int j = i+1;
        float t = pos - i;
        RecordSampleType& si = buffer[i&record_buffer_size_mask_];
        RecordSampleType& sj = buffer[j&record_buffer_size_mask_];
        RecordSampleType grn = vessl::math::lerp(si, sj, t) * env;
        vessl::sample::spatialize(grn.value(), grain.pan, &samp);
        accum += samp;
      }
      
      grain.ramp += grain.speed;
      
      // swap with last active grain when finished
      if (grain.ramp >= grain.size)
      {
        grains_[i] = grains_[--active_grain_count_];
      }
    }
    
    return accum;
  }
  
  VESSL_INLINE void generate(vessl::array<SampleType> out)
  {
    if (const float grain_spacing = params_.grain_rate.value.samples; grain_spacing > 0)
    {
      float trigger_delay = grain_spacing - grain_rate_phasor_;
      grain_rate_phasor_ += out.size();
    
      while (grain_rate_phasor_ >= grain_spacing)
      {
        grain_rate_phasor_ -= grain_spacing;
        trigger(trigger_delay);
        trigger_delay += grain_spacing;
      }
    }
    
    out.fill(SampleType(0));
    
    SampleType samp;
    for (int g = active_grain_count_ - 1; g >= 0; g--)
    {
      Grain& grain = grains_[g];
      const float grain_pos = grain.start + grain.ramp;
      // block copy the number of samples we'll need for this grain from our buffer into our scratch space
      vessl::size_t block_size = out.size() * 2;
      vessl::size_t scratch_start = static_cast<size_t>(grain_pos);
      vessl::size_t scratch_end = scratch_start + block_size;

      if (block_size > 0)
      {
        scratch_start &= record_buffer_size_mask_;
        scratch_end &= record_buffer_size_mask_;
        if (scratch_start < scratch_end)
        {
          vessl::array scratch(scratch_buffer_, block_size);
          vessl::array block(record_buffer_.data() + scratch_start, block_size);
          block.copy_to(scratch);
        }
        else
        {
          vessl::size_t to_end = record_buffer_.size() - scratch_start;
          
          vessl::array block_a(record_buffer_.data() + scratch_start, to_end);
          vessl::array scratch_a(scratch_buffer_, to_end);
          block_a.copy_to(scratch_a);
          
          vessl::array block_b(record_buffer_.data(), block_size - to_end);
          vessl::array scratch_b(scratch_buffer_ + to_end, block_size - to_end);
          block_b.copy_to(scratch_b);
        }
      
        float scratch_pos = grain_pos - vessl::math::floor(grain_pos);
        bool grain_done = false;
        for (int i = 0; i < out.size() && !grain_done; i++)
        {
          SampleType& gro = out[i];
          if (grain.ramp >= 0)
          {
            float env = grain.ramp < grain.decay_start 
                      ? grain.ramp * grain.attack_mult 
                      : (grain.size - grain.ramp) * grain.decay_mult;
            int x = static_cast<int>(scratch_pos);
            int y = x+1;
            float t = scratch_pos - x;
            RecordSampleType& si = scratch_buffer_[x];
            RecordSampleType& sj = scratch_buffer_[y];
            RecordSampleType grn = vessl::math::lerp(si, sj, t) * env;
            vessl::sample::spatialize(grn.value(), grain.pan, &samp);
            gro += samp;
          }
          
          scratch_pos += grain.speed;
          grain.ramp += grain.speed;
          
          grain_done = grain.ramp >= grain.size;
        }
        
        // swap with last active grain when finished
        if (grain_done)
        {
          grains_[g] = grains_[--active_grain_count_];
        }
      }
    }
  }
  
  // buffer_size must be a power of two!
  static Granulator* create(vessl::size_t buffer_size, vessl::size_t block_size)
  {
    RecordSampleType* scratch = new RecordSampleType[block_size*2];
    RecordSampleType* buffer = new RecordSampleType[buffer_size];
    Granulator* granulator = new Granulator(buffer, buffer_size);
    granulator->scratch_buffer_ = scratch;
    return granulator;
  }
  
  static void destroy(Granulator* granulator)
  {
    delete[] granulator->scratch_buffer_;
    delete[] granulator->record_buffer_.data();
    delete granulator;
  }

protected:
  VESSL_INLINE Parameter element_at(vessl::size_t index) const override
  {
    switch (index)
    {
      case 0: return grain_duration();
      case 1: return grain_speed();
      case 2: return grain_offset();
      case 3: return grain_rate();
      case 4: return grain_pan();
      default: return Parameter::none();
    }
  }
  
private:
  Granulator(RecordSampleType* buffer, size_t buffer_size)
    : record_buffer_(buffer, buffer_size)
    , record_buffer_size_mask_(buffer_size - 1)
    , scratch_buffer_(nullptr)
  {

  }

  struct 
  {
    vessl::duration_p grain_duration;
    vessl::analog_p   grain_speed;
    vessl::duration_p grain_offset;
    vessl::duration_p grain_rate;
    vessl::analog_p   grain_pan;
  } params_;
  
  struct Grain
  {
    float start;
    float size;
    float speed;
    float pan;
    float attack_mult;
    float decay_start;
    float decay_mult;
    
    float ramp;
  };
  
  unsigned grain_rate_phasor_ = 0;
  unsigned active_grain_count_ = 0;
  unsigned grain_triggered_ = 0;
  Grain grains_[MaxGrains];
  
  vessl::sample::delay_line<RecordSampleType> record_buffer_;
  vessl::size_t record_buffer_size_mask_;
  
  RecordSampleType* scratch_buffer_;
};