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
  using Freeze   = vessl::processors::freeze<T>;
  using BitCrush = vessl::processors::bitcrush<T, 24>;
  using Ramp     = vessl::generators::ramp<analog_t>;
  using Parameter = vessl::parameter;

public:
  Cutup(analog_t sample_rate, T* freeze_buffer, size_t freeze_buffer_size)
  : clockable(sample_rate, 16, sample_rate / 4, 120.f)
  , freeze_(vessl::array(freeze_buffer, freeze_buffer_size), sample_rate)
  , crush_(sample_rate, sample_rate)
  , slice_dur_(0)
  , freeze_enabled_(freeze_.enabled())
  , slice_amp_(1.f)
  , slice_damp_(1.f)
  , rate_ramp_(sample_rate)
  , slice_ramp_(sample_rate)
  , crush_ramp_(sample_rate)
  {
    memset(freeze_buffer, 0, sizeof(T)*freeze_buffer_size);
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
    binary_t frozen = freeze_enabled_.read<binary_t>();
    if (enable && !frozen)
    {
      freeze_enabled_ = true;
      start_slice();
    }
    else if (!enable && frozen)
    {
      freeze_enabled_ = false;
      slice_amp_ = 1.f;
      set_crush(0);
    }
    else if (frozen)
    {
      analog_t dur = slice_ramp_.generate();
      freeze_.position() = slice_start_ - dur;
      freeze_.duration() = duration_t(dur);
      freeze_.rate() = rate_ramp_.generate<vessl::math::easing::expo::out>();
      set_crush(crush_ramp_.generate());
      if (freeze_.looped())
      {
        slice_amp_ *= slice_damp_;
      }
      if (--slice_dur_.samples <= 0)
      {
        start_slice();
      }
    }
    sample_t out;
    in >> freeze_ >> crush_ >> out;
    return out * slice_amp_;
  }

protected:
  void set_crush(analog_t amt)
  {
    float sr   = clockable::sample_rate_;
    float bits = amt > 0.001f ? (16.f - amt*12.0f) : 24;
    float rate = amt > 0.001f ? sr * 0.25f + amt*(100 - sr * 0.25f) : sr;
    crush_.depth() = bits;
    crush_.rate() = rate;
  }

  void start_slice()
  {
    // 2 ticks is half a beat in our reckoning
    int slice_ticks = 2;
    slice_dur_      = duration_t(clockable::tempo_.samples * slice_ticks);
    // we choose one of the previous eight slice_dur_ sections to slice (i.e. a section from the last 4/4 measure)
    int slice_idx = vessl::math::random::range(1, 9);
    int slice_clf = 1 << vessl::math::random::range(params_.slice_min.value, params_.slice_max.value);
    int slice_clt = 1 << vessl::math::random::range(params_.slice_min.value, params_.slice_max.value);
    int rate_min  = vessl::math::max(params_.rate_min.value, 1);
    int slice_rtf = vessl::math::random::range(rate_min, params_.rate_max.value+1);
    int slice_rtt = vessl::math::random::range(params_.rate_min.value, params_.rate_max.value+1);
    int slice_dir = vessl::math::random::range(0, 2) ? 1 : -1;
    analog_t crush_f = vessl::math::random::range(params_.crush_min.value, params_.crush_max.value);
    analog_t crush_t = vessl::math::random::range(params_.crush_min.value, params_.crush_max.value); 

    slice_start_ = slice_idx * slice_dur_.samples;
    slice_amp_   = 1.f;
    slice_damp_  = 1.0f - vessl::math::random::range(0.f, params_.damp_max.value);

    analog_t slice_dur_f = clockable::tempo_.samples * 0.25f * slice_clf;
    analog_t slice_dur_t = clockable::tempo_.samples * 0.25f * slice_clt;
    
    freeze_.position()   = slice_start_ - slice_dur_f;
    freeze_.duration()   = duration_t(slice_dur_f);
    freeze_.rate()       = static_cast<float>(slice_rtf) * slice_dir;
    freeze_.reset();
    
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
  // processors
  Freeze freeze_;
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

  duration_t slice_dur_;
  analog_t   slice_start_;
  analog_t   slice_amp_;
  analog_t   slice_damp_;

  Parameter freeze_enabled_; // cache freeze's enabled parameter because we need to access it so much.
  Ramp rate_ramp_;
  Ramp slice_ramp_;
  Ramp crush_ramp_;
};
