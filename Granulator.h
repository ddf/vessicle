#pragma once

#include "vessl/vessl.h"

template<typename T, unsigned ChannelCount, unsigned MaxGrains>
class Granulator : public vessl::unit_processor<vessl::sample::frame<T, ChannelCount>>
                 , public vessl::generator<vessl::sample::frame<T, ChannelCount>>
                 , protected vessl::plist<7>
{
public:
  using Parameter = vessl::parameter;
  using SampleType = vessl::sample::frame<T, ChannelCount>;
  using RecordSampleType = typename vessl::sample::type<T>::mono;
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
  
  GrainEnvelope envelope;
  
  // starts a new grain
  VESSL_INLINE void trigger(float sample_delay = 0)
  {
    if (active_grain_count_ < MaxGrains)
    {
      Grain& grain = grains_[active_grain_count_++];
      grain.speed = vessl::math::constrain(params_.grain_speed.value, 0.5f, 2.f);
      grain.dir = params_.grain_reverse.value ? -1.f : 1.f;
      // start is always relative to the end of buffer_a (i.e. the beginning of buffer_b)
      // so that no matter the grain size or offset, we will always read from a valid portion
      // of the double-sized buffer shared by the two delay lines.
      if (params_.grain_reverse.value)
      {
        // start at the end of the grain and play in reverse
        grain.size = params_.grain_duration.value.samples * grain.speed * -1.f;
        grain.start = record_buffer_a_.size() 
                    + record_buffer_a_.get_write_index()
                    - params_.grain_offset.value.samples;
      }
      else
      {
        grain.size = params_.grain_duration.value.samples * grain.speed;
        grain.start = record_buffer_a_.size() 
                    + record_buffer_a_.get_write_index()
                    - grain.size 
                    - params_.grain_offset.value.samples;
      }
      grain.pan = vessl::math::constrain(params_.grain_pan.value, -1.f, 1.f);
      grain.vol = params_.grain_volume.value;
      grain.ramp = sample_delay > 0 ? -(static_cast<float>(vessl::phase_360) / sample_delay) : 0;
      grain.ramp_step = vessl::phase_360 / params_.grain_duration.value.samples;
      
      grain_triggered_ = true;
    }
  }
    
  // should only be called immediately after process/generate
  [[nodiscard]] bool started_grain() const { return grain_triggered_; }
  
  [[nodiscard]] int active_grain_count() const { return active_grain_count_; }

  [[nodiscard]] VESSL_INLINE SampleType process(const SampleType &in) override
  {
    grain_triggered_ = false;
    RecordSampleType rin = in.to_mono();
    record_buffer_a_.write(rin);
    record_buffer_b_.write(rin);
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
      RecordSampleType rs = s.to_mono();
      record_buffer_a_.write(rs);
      record_buffer_b_.write(rs);
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
    RecordSampleType* buffer = record_buffer_a_.data();
    constexpr vessl::analog_t to_analog = 1.0f / vessl::phase_360;
    for (int i = active_grain_count_ - 1; i >= 0; i--)
    {
      Grain& grain = grains_[i];
      if (grain.ramp >= 0)
      {
        float gt  = grain.ramp*to_analog;
        float pos = grain.start + grain.size * gt;
        T env = envelope.evaluate(static_cast<vessl::phase_t>(grain.ramp)) * grain.vol;
        RecordSampleType grn = vessl::sample::readf<InterpType>(buffer, pos) * env;
        vessl::sample::spatialize(grn.value(), grain.pan, &samp);
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
        // because we double buffer the record buffer 
        // and set grain.start relative to the middle of that buffer,
        // we can always read a continguous block!
        vessl::array scratch(scratch_buffer_, scratch_write_size);
        vessl::array block(record_buffer_a_.data() + scratch_start, scratch_write_size);
        block.copy_to(scratch);
      
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
            RecordSampleType grn = vessl::sample::readf<InterpType>(scratch_buffer_, scratch_pos) * env;
            vessl::sample::spatialize(grn.value(), grain.pan, &samp);
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
    RecordSampleType* scratch = new RecordSampleType[block_size*2];
    RecordSampleType* buffer = new RecordSampleType[buffer_size*2];
    Granulator* granulator = new Granulator(buffer, buffer_size);
    granulator->scratch_buffer_ = scratch;
    return granulator;
  }
  
  static void destroy(Granulator* granulator)
  {
    delete[] granulator->scratch_buffer_;
    delete[] granulator->record_buffer_a_.data();
    // don't need to delete record buffer b data because a+b were allocated as a contiguous block.
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
      default: return Parameter::none();
    }
  }
  
private:
  Granulator(RecordSampleType* buffer, size_t buffer_size)
    : record_buffer_a_(buffer, buffer_size)
    , record_buffer_b_(buffer + buffer_size, buffer_size)
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
    vessl::analog_p   grain_volume;
    vessl::binary_p   grain_reverse;
  } params_;
  
  struct Grain
  {
    vessl::digital_t ramp; // phase, but allowing for negative
    vessl::digital_t ramp_step;
    
    float start;
    float size;
    float speed;
    float pan;
    float vol;
    float dir; // +1 forward, -1 backward
  };

  Grain grains_[MaxGrains];
  
  // we allocate twice as much buffer as we need and split it between two delay lines.
  // this allows us to always block copy a contiguous block of sample data
  // and also allows us to use vessl::sample::read_interpolated with our buffer data
  // without needing to do anything special with reads between the end and beginning of the buffer.
  RecordSampleType* scratch_buffer_;
  vessl::sample::delay_line<RecordSampleType> record_buffer_a_;
  vessl::sample::delay_line<RecordSampleType> record_buffer_b_;
  
  unsigned grain_rate_phasor_ = 0;
  unsigned active_grain_count_ = 0;
  unsigned grain_triggered_ = 0;
};