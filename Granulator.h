#pragma once

#include "vessl/vessl.h"

template<typename T, unsigned ChannelCount, unsigned MaxGrains>
class Granulator : public vessl::unit_processor<vessl::sample::frame<T, ChannelCount>>
                 , public vessl::generator<vessl::sample::frame<T, ChannelCount>>
                 , protected vessl::plist<8>
{
public:
  using Parameter = vessl::parameter;
  using SampleType = vessl::sample::frame<T, ChannelCount>;
  using InterpType = vessl::sample::interpolation::linear;
  using GrainEnvelope = vessl::sample::waves::unipolar::triangle<T>;
  
  [[nodiscard]] const parameter_list& parameters() const override { return *this; }
  
  // length of a grain in seconds (duration_t)
  [[nodiscard]] Parameter duration() const { return params_.grain_duration("g.dur", 'd'); }
  // playback speed of grain, clamped to [-2,2]
  [[nodiscard]] Parameter speed() const { return params_.grain_speed("g.spd", 's'); }
  // how far back in time the grain start should be moved (duration_t)
  [[nodiscard]] Parameter offset() const { return params_.grain_offset("g.off", 'o'); }
  // how often grains are started (duration_t)
  [[nodiscard]] Parameter rate() const { return params_.grain_rate("g.rate", 'r'); }
  // how a grain is panned in the multi-channel field [-1,1]
  [[nodiscard]] Parameter pan() const { return params_.grain_pan("g.pan", 'p'); }
  // linear scale on amplitude of grain
  [[nodiscard]] Parameter volume() const { return params_.grain_volume("g.vol", 'v'); }
  // grains will play in reverse if true
  [[nodiscard]] Parameter reverse() const { return params_.grain_reverse("g.rev", 'f'); }
  // the maximum number of grains that can be simultaneously active.
  [[nodiscard]] Parameter max_active() const { return params_.max_grains("g.max", 'm'); }
  
  GrainEnvelope envelope;
  
  // starts a new grain
  VESSL_INLINE void trigger(float sample_delay = 0)
  {
    unsigned max_active = vessl::math::min(params_.max_grains.value, MaxGrains);
    if (active_grain_count_ < max_active)
    {
      Grain& grain = grains_[active_grain_count_++];
      grain.speed = vessl::math::constrain(params_.grain_speed.value, 0.5f, 2.f);
      grain.dir = params_.grain_reverse.value ? -1.f : 1.f;
      // start is always relative to the end of the buffer
      // so that no matter the grain size or offset,
      // we will always read with a positive index.
      if (params_.grain_reverse.value)
      {
        // start at the end of the grain and play in reverse
        grain.size = params_.grain_duration.value.samples * grain.speed * -1.f;
        grain.start = record_buffer_.size() 
                    + record_buffer_.get_write_index()
                    - params_.grain_offset.value.samples;
      }
      else
      {
        grain.size = params_.grain_duration.value.samples * grain.speed;
        grain.start = record_buffer_.size() 
                    + record_buffer_.get_write_index()
                    - grain.size 
                    - params_.grain_offset.value.samples;
      }
      float bal = vessl::math::constrain(params_.grain_pan.value, -1.f, 1.f);
      vessl::sample::make_spatializer(bal, &grain.mix);
      grain.vol = params_.grain_volume.value;
      grain.ramp = sample_delay > 0 ? -(static_cast<float>(vessl::phase_360) / sample_delay) : 0;
      grain.ramp_step = vessl::phase_360 / params_.grain_duration.value.samples;
      
      grain_triggered_ = true;
    }
  }
    
  // should only be called immediately after process/generate
  [[nodiscard]] VESSL_INLINE bool started_grain() const { return grain_triggered_; }
  
  [[nodiscard]] VESSL_INLINE int active_grain_count() const { return active_grain_count_; }

  [[nodiscard]] VESSL_INLINE SampleType process(const SampleType &in) override
  {
    record_buffer_.write(in);
    return generate();
  }
  
  VESSL_INLINE void process(const vessl::array<SampleType>& in, vessl::array<SampleType> out)
  {
    auto rin = in.make_reader();
    //auto wout = out.make_writer();
    while (rin)
    {
      SampleType s = rin.read();
      record_buffer_.write(s);
      //wout << generate();
    }
    
    generate(out);
  }
  
  VESSL_INLINE void overdub(vessl::source<SampleType>& in, vessl::analog_t scale, size_t sample_delay = 0)
  {
    size_t write_offset = sample_delay 
                        + record_buffer_.size()
                        - params_.grain_offset.value.samples
                        - params_.grain_duration.value.samples;
    while (!in.is_empty())
    {
      SampleType rs = in.read() * scale;
      SampleType od = record_buffer_.overdub(rs, write_offset);
      ++write_offset;
    }
  }
  
  [[nodiscard]] VESSL_INLINE SampleType generate() override
  {
    grain_triggered_ = false;
    grain_rate_phasor_++;
    
    if (grain_rate_phasor_ >= params_.grain_rate.value.samples)
    {
      trigger();
      grain_rate_phasor_ = 0;
    }
    
    SampleType accum = SampleType(0);
    SampleType samp;
    SampleType* buffer = record_buffer_.data();
    constexpr vessl::analog_t to_analog = 1.0f / vessl::phase_360;
    for (int i = active_grain_count_ - 1; i >= 0; i--)
    {
      Grain& grain = grains_[i];
      if (grain.ramp >= 0)
      {
        float gt  = grain.ramp*to_analog;
        float pos = grain.start + grain.size * gt;
        T env = envelope.evaluate(static_cast<vessl::phase_t>(grain.ramp)) * grain.vol;
        SampleType grn = vessl::sample::readf<InterpType>(buffer, pos) * env;
        grain.mix.spatialize(grn, &samp);
        accum += samp;
      }
      
      grain.ramp += grain.ramp_step;
      
      // swap with last active grain when finished
      if (grain.ramp >= vessl::phase_360)
      {
        grains_[i] = grains_[--active_grain_count_];
      }
    }
    
    return accum;
  }
  
  VESSL_INLINE void generate(vessl::array<SampleType> out)
  {
    grain_triggered_ = false;
    
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
    constexpr vessl::analog_t to_analog = 1.0f / static_cast<float>(vessl::phase_360);
    SampleType* buffer = record_buffer_.data();
    for (int g = active_grain_count_ - 1; g >= 0; g--)
    {
      Grain& grain = grains_[g];
      const float grain_pos = grain.start + grain.size * (grain.ramp * to_analog);
      // block copy the maximum number of samples we'll need for this grain from our buffer into our scratch space
      vessl::size_t scratch_write_size = out.size() * 2;
      vessl::size_t scratch_start = grain.dir > 0 
        ? static_cast<size_t>(grain_pos) 
        : static_cast<size_t>(grain_pos + 2) - scratch_write_size;
      vessl::size_t scratch_end = scratch_start + scratch_write_size;

      if (scratch_write_size > 0)
      {
        scratch_start &= record_buffer_size_mask_;
        scratch_end &= record_buffer_size_mask_;
        
        if (scratch_start < scratch_end)
        {
          vessl::array scratch(scratch_buffer_, scratch_write_size);
          vessl::array block(buffer + scratch_start, scratch_write_size);
          block.copy_to(scratch);
        }
        else
        {
          vessl::size_t block_a_size = record_buffer_size_mask_ + 1 - scratch_start;
          vessl::array block_a(buffer + scratch_start, block_a_size);
          vessl::array scratch_a(scratch_buffer_, block_a_size);
          vessl::array block_b(buffer, scratch_write_size - block_a_size);
          vessl::array scratch_b(scratch_buffer_ + block_a_size, scratch_write_size - block_a_size);
          block_a.copy_to(scratch_a);
          block_b.copy_to(scratch_b);
        }
      
        float scratch_pos = grain.dir > 0 
          ? grain_pos - vessl::math::floor(grain_pos)
          : static_cast<float>(scratch_write_size - 2) + (grain_pos - vessl::math::floor(grain_pos));
        float scratch_step = grain.speed*grain.dir;
        bool grain_done = false;
        for (int i = 0; i < out.size() && !grain_done; i++)
        {
          SampleType& gro = out[i];
          if (grain.ramp >= 0)
          {
            T env = envelope.evaluate(static_cast<vessl::phase_t>(grain.ramp)) * grain.vol;
            SampleType grn = vessl::sample::readf<InterpType>(scratch_buffer_, scratch_pos) * env;
            grain.mix.spatialize(grn, &samp);
            gro += samp;
          }
          
          scratch_pos += scratch_step;
          grain.ramp += grain.ramp_step;
          
          grain_done = grain.ramp >= vessl::phase_360;
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
    SampleType* scratch = new SampleType[block_size*2];
    SampleType* buffer = new SampleType[buffer_size];
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
      case 0: return duration();
      case 1: return speed();
      case 2: return offset();
      case 3: return rate();
      case 4: return pan();
      case 5: return volume();
      case 6: return reverse();
      case 7: return max_active();
      
      default: return Parameter::none();
    }
  }
  
private:
  Granulator(SampleType* buffer, size_t buffer_size)
    : record_buffer_(buffer, buffer_size)
    , scratch_buffer_(nullptr)
    , record_buffer_size_mask_(buffer_size - 1)
  {

  }

  struct 
  {
    vessl::duration_p grain_duration;
    vessl::analog_p   grain_speed;
    vessl::duration_p grain_offset;
    vessl::duration_p grain_rate;
    vessl::analog_p   grain_pan;
    vessl::analog_p   grain_volume;
    vessl::binary_p   grain_reverse;
    vessl::param<unsigned> max_grains;
  } params_;
  
  struct Grain
  {
    vessl::digital_t ramp; // phase, but allowing for negative
    vessl::digital_t ramp_step;
    
    float start;
    float size;
    float speed;
    float vol;
    float dir; // +1 forward, -1 backward
    SampleType mix;
  };

  Grain grains_[MaxGrains];
  
  vessl::sample::delay_line<SampleType> record_buffer_;
  // to improve performance, we block copy audio from the record buffer 
  // to a smaller scratch buffer before generating each active grain.
  SampleType* scratch_buffer_;
  
  unsigned record_buffer_size_mask_;
  unsigned grain_rate_phasor_ = 0;
  unsigned active_grain_count_ = 0;
  unsigned grain_triggered_ = 0;
};