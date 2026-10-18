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

template<typename T, size_t BufferSize>
class Cutup : public vessl::unit_processor<T>
            , public vessl::time::clockable
            , public vessl::plist<1>
{
  using binary_t = vessl::binary_t;
  using analog_t = vessl::analog_t;
  using duration_t = vessl::time::duration;
  using size_t   = vessl::size_t;
  using sample_t = T;
  using Freeze   = vessl::processors::freeze<T>;
  using Parameter = vessl::parameter;

public:
  Cutup(analog_t sample_rate)
  : clockable(sample_rate, 16, sample_rate / 4, 120.f)
  , slice_dur_(0)
  , freeze_(vessl::array(freeze_buffer_, BufferSize), sample_rate)
  , freeze_enabled_(freeze_.enabled())
  {
    memset(freeze_buffer_, 0, sizeof(T)*BufferSize);
  }

  Parameter enabled() const { return params_.enabled("enabled", 'e'); }

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
    }
    else if (frozen)
    {
      if (--slice_dur_.samples <= 0)
      {
        start_slice();
      }
    }
    return freeze_.process(in);
  }

protected:
  void start_slice()
  {
    int slice_ticks = 2 * vessl::math::random::range(1, 3);
    slice_dur_      = duration_t(clockable::tempo_.samples * slice_ticks);
    // we one of the previous four beats to slice
    int slice_idx = vessl::math::random::range(1, 5);
    int slice_clk = vessl::math::random::range(1, 5);
    int slice_rte = vessl::math::random::range(1, 9);
    int slice_dir = vessl::math::random::range(0, 2) ? 1 : -1;
    analog_t slice_start = slice_idx * slice_dur_.samples;
    analog_t slice_end   = slice_start - clockable::tempo_.samples * slice_clk;
    freeze_.position()   = slice_end;
    freeze_.duration()   = duration_t(slice_start - slice_end);
    freeze_.rate()       = static_cast<float>(slice_rte) * slice_dir;
  }

  [[nodiscard]] Parameter element_at(size_t index) const override
  {
    switch(index)
    {
      case 0: return enabled();
    }

    return Parameter::none();
  }

private:
  sample_t freeze_buffer_[BufferSize];

  struct
  {
    vessl::binary_p enabled;
  } params_;

  duration_t slice_dur_;

  Freeze freeze_;
  Parameter freeze_enabled_; // cache freeze's enabled parameter because we need to access it so much.
};
