/**

AUTHOR:
    (c) 2026 Damien Quartz

LICENSE:
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


DESCRIPTION:

*/

#pragma once

#include "vessl/vessl.h"

template<typename T>
class Cutup : public vessl::unit_processor<T>
            , public vessl::time::clockable
            , public vessl::plist<5>
{
  using binary_t = vessl::binary_t;
  using analog_t = vessl::analog_t;
  using duration_t = vessl::time::duration;
  using size_t   = vessl::size_t;
  using sample_t = T;
  using Buffer   = vessl::sample::ring_buffer<T>;
  using BitCrush = vessl::processors::bitcrush<T, 24>;
  using Ramp     = vessl::generators::ramp<analog_t>;
  using Parameter = vessl::parameter;

  static constexpr int slice_count    = 16;
  static constexpr int slice_tick_len = 2;

public:
  Cutup(analog_t sample_rate, T* buffer, size_t buffer_size)
  : clockable(sample_rate, 16, sample_rate / 4, 120.f)
  , buffer_(buffer, buffer_size)
  , crush_(sample_rate, sample_rate)
  , slices_head_(slice_count - 1)
  , slices_tick_(slice_tick_len)
  , slice_enabled_(false)
  , slice_dur_(0)
  , slice_amp_(1.f)
  , slice_damp_(1.f)
  , rate_ramp_(sample_rate)
  , slice_ramp_(sample_rate)
  , crush_ramp_(sample_rate)
  {
    memset(buffer, 0, sizeof(T)*buffer_size);
    params_.rate_min.value = 1;
    params_.rate_max.value = 1;
    params_.slice_min.value = 0;
    params_.slice_max.value = 6;
    set_crush(0);
  }

  Parameter enabled() const { return params_.enabled("enabled", 'e'); }
  Parameter slice_min() const { return params_.slice_min("slice min", 's'); }
  Parameter slice_max() const { return params_.slice_max("slice max", 'S'); }
  // range for randomly chosen rate per slice (inclusive)
  Parameter rate_min() const { return params_.rate_min("rate min", 'r'); }
  Parameter rate_max() const { return params_.rate_max("rate max", 'R'); }
  // range for randomly chosen bit crush per slice [0,1] where 0 is no bit crushing.
  Parameter crush_min() const { return params_.crush_min("crush min", 'c'); }
  Parameter crush_max() const { return params_.crush_max("crush max", 'C'); }
  // range for how quickly single slice is damped over time [0,1]
  Parameter damp_max() const { return params_.damp_max("damp max", 'D'); }

  const vessl::parameter_list& parameters() const override { return *this; }

  sample_t process(const sample_t& in) override
  {
    clockable::tick();
    binary_t enable = params_.enabled.value;
    if (enable && !slice_enabled_)
    {
      slice_enabled_ = true;
      start_slice();
    }
    else if (!enable && slice_enabled_)
    {
      slice_enabled_ = false;
      slice_amp_ = 1.f;
      set_crush(0);
    }
    else if (slice_enabled_)
    {
      analog_t dur = slice_ramp_.generate();
      slice_end_ = slice_start_ + dur;
      slice_rate_ = rate_ramp_.generate<vessl::math::easing::expo::out>();
      set_crush(crush_ramp_.generate());
      if (--slice_dur_.samples <= 0)
      {
        start_slice();
      }
    }
    sample_t out;
    process_slice(in) >> crush_ >> out;
    return out * slice_amp_;
  }

protected:
  VESSL_INLINE void tock(period_t sample_delay) override
  {
    if (!slice_enabled_)
    {
      int len = vessl::math::max(static_cast<int>(tempo_.samples), 4);
      if (slices_tick_ == slice_tick_len)
      {
        slices_head_ = (slices_head_+1)%slice_count;
        slice& slice = slices_[slices_head_];

        slice.idx = (buffer_.get_write_index() + sample_delay + buffer_.size() - len) % buffer_.size();
        slice.len = len;
        slices_tick_ = 0;
      }
      else
      {
        slices_[slices_head_].len += len;
        slices_tick_++;
      }
    }
  }

  VESSL_INLINE sample_t process_slice(const sample_t& in)
  {
    sample_t out = in;
    if (slice_enabled_)
    {
      T* data = buffer_.data();
      size_t sz = buffer_.size();
      size_t x = static_cast<size_t>(slice_pos_);
      size_t y = x+1;
      const sample_t& a = data[(x+sz)%sz];
      const sample_t& b = data[(y+sz)%sz]; 
      analog_t        t = slice_pos_ - x;
      out = a + (b-a)*t;
      
      slice_pos_ += slice_rate_;
      const analog_t len = slice_end_ - slice_start_;
      if (slice_pos_ < slice_start_)
      {
        slice_pos_ += len;
        slice_amp_ *= slice_damp_;
      }
      else if (slice_pos_ > slice_end_)
      {
        slice_pos_ -= len;
        slice_amp_ *= slice_damp_;
      }
    }
    else
    {
      buffer_.write(in);
    }
    return out;
  }

  VESSL_INLINE void set_crush(analog_t amt)
  {
    float sr   = clockable::sample_rate_;
    float bits = amt > 0.001f ? (16.f - amt*12.0f) : 24;
    float rate = amt > 0.001f ? sr * 0.25f + amt*(100 - sr * 0.25f) : sr;
    crush_.depth() = bits;
    crush_.rate() = rate;
  }

  VESSL_INLINE void start_slice()
  {
    // 2 ticks is half a beat in our reckoning
    static constexpr int st[] = { 4, 6, 8, 12 };
    static constexpr int stc = sizeof(st) / sizeof(int);
    int slice_ticks = st[vessl::math::random::range(0, stc)];
    // for how long in samples this slice should be played back
    slice_dur_      = duration_t(clockable::tempo_.samples * slice_ticks);
    // we choose one of the previous eight slice_dur_ sections to slice (i.e. a section from the last 4/4 measure)
    int slice_idx = vessl::math::random::range(0, slice_count);
    int slice_clf = 1 << vessl::math::random::range(params_.slice_min.value, params_.slice_max.value);
    int slice_clt = 1 << vessl::math::random::range(params_.slice_min.value, params_.slice_max.value);
    int rate_min  = vessl::math::max(params_.rate_min.value, 1);
    int slice_rtf = vessl::math::random::range(rate_min, params_.rate_max.value+1);
    int slice_rtt = vessl::math::random::range(params_.rate_min.value, params_.rate_max.value+1);
    int slice_dir = vessl::math::random::range(0, 2) ? 1 : -1;
    analog_t crush_f = vessl::math::random::range(params_.crush_min.value, params_.crush_max.value);
    analog_t crush_t = vessl::math::random::range(params_.crush_min.value, params_.crush_max.value); 

    const slice& slice = slices_[slice_idx];
    slice_start_ = slice.idx;
    slice_amp_   = 1.f;
    slice_damp_  = 1.0f - vessl::math::random::range(0.f, params_.damp_max.value);

    analog_t slice_dur_f = vessl::math::max(2.f, clockable::tempo_.samples * 0.25f * slice_clf);
    analog_t slice_dur_t = vessl::math::max(2.f, clockable::tempo_.samples * 0.25f * slice_clt);
    
    slice_end_           = slice_start_ + slice_dur_f;
    slice_rate_          = static_cast<float>(slice_rtf) * slice_dir;
    slice_pos_           = slice_dir > 0 ? slice_start_ : slice_end_;
    
    rate_ramp_.from()     = static_cast<float>(slice_rtf) * slice_dir;
    rate_ramp_.to()       = static_cast<float>(slice_rtt) * slice_dir;
    rate_ramp_.duration() = slice_dur_.to_seconds(clockable::sample_rate_);
    rate_ramp_.trigger();

    slice_ramp_.from()     = slice_dur_f;
    slice_ramp_.to()       = slice_dur_t;
    slice_ramp_.duration() = slice_dur_.to_seconds(clockable::sample_rate_);
    slice_ramp_.trigger();

    crush_ramp_.from()     = crush_f;
    crush_ramp_.to()       = crush_t;
    crush_ramp_.duration() = slice_dur_.to_seconds(clockable::sample_rate_);
    crush_ramp_.trigger();

    set_crush(crush_f);
  }

  [[nodiscard]] Parameter element_at(size_t index) const override
  {
    switch(index)
    {
      case 0: return enabled();
      case 1: return rate_min();
      case 2: return rate_max();
      case 3: return slice_min();
      case 4: return slice_max();
      case 5: return crush_min();
      case 6: return crush_max();
      case 7: return damp_max();
    }

    return Parameter::none();
  }

private:
  Buffer buffer_;
  // processors
  BitCrush crush_;

  // paramers and state
  struct
  {
    vessl::binary_p enabled;
    vessl::param<int> rate_min;
    vessl::param<int> rate_max;
    vessl::param<int> slice_min;
    vessl::param<int> slice_max;
    vessl::analog_p crush_min;
    vessl::analog_p crush_max;
    vessl::analog_p damp_max;
  } params_;

  struct slice
  {
    // index in buffer_ where this slice begins
    size_t idx = 0;
    // length in samples of the slice 
    // (based on value of tempo_duration.samples at the time of recording)
    size_t len = 2;
  };
  
  slice slices_[slice_count];
  int slices_head_, slices_tick_;

  binary_t   slice_enabled_;
  analog_t   slice_start_;
  analog_t   slice_end_;
  analog_t   slice_pos_;
  duration_t slice_dur_;
  analog_t   slice_rate_;
  analog_t   slice_amp_;
  analog_t   slice_damp_;

  Ramp rate_ramp_;
  Ramp slice_ramp_;
  Ramp crush_ramp_;
};
