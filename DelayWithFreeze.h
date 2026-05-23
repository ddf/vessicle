#pragma once

#include "vessl/vessl.h"

template<typename T>
class DelayWithFreeze : public vessl::unit_processor<T>, protected vessl::plist<5>
{
public:
  using sample_t = T;
  using param = vessl::parameter;
  using binary_t = vessl::binary_t;
  using analog_t = vessl::analog_t;
  using array = vessl::array<T>;
  
  DelayWithFreeze(array buffer, analog_t sample_rate, analog_t delay_in_seconds = 0, analog_t feedback = 0)
    : vessl::unit_processor<T>()
    , fader_(0.95f, 0)
    , delay_(buffer, sample_rate, delay_in_seconds, feedback)
    , freeze_(buffer, sample_rate)
  {
  }

  [[nodiscard]] const vessl::parameter_list& parameters() const override { return *this; }

  [[nodiscard]] param time() const { return delay_.time(); }
  [[nodiscard]] param feedback() const { return delay_.feedback(); }
  [[nodiscard]] param freeze_enabled() const
  {
    return params_.frozen("freeze enabled", 'e');
  }
  [[nodiscard]] param freeze_position() const { return freeze_.position(); }
  [[nodiscard]] param freeze_duration() const { return freeze_.duration(); }

  sample_t process(const T& in) override
  {
    binary_t frozen = params_.frozen.value;
    analog_t fade = fader_ = (frozen ? 1.0f : 0.0f);
    sample_t s1 = frozen ? in : delay_.process(in);
    if (!frozen)
    {
      freeze_.buffer().set_write_index(delay_.buffer().get_write_index());
    } 
    T s2 = fade > 0 ? freeze_.generate() : 0.f;
    return vessl::sample::crossfade(s1, s2, fade);
  }

  template<vessl::time::mode TimeMode = vessl::time::mode::slew>
  void process(array input, array output)
  {
    if (params_.frozen.value)
    {
      freeze_.buffer().set_write_index(delay_.buffer().get_write_index());
      if (fader_.value < 0.999f)
      {
        auto r = input.make_reader();
        auto w = output.make_writer();
        while (r)
        {
          w << process(r.read());
        }
      }
      else
      {
        freeze_.template generate<TimeMode>(output);
      }
    }
    else
    {
      if (fader_.value > 0.001f)
      {
        auto r = input.make_reader();
        auto w = output.make_writer();
        while (r)
        {
          w << process(r.read());
        }
      }
      else
      {
        delay_.template process<TimeMode>(input, output);
      }
    }
  }
  
protected:
  [[nodiscard]] param element_at(vessl::size_t index) const override
  {
    param p[num] = { time(), feedback(), freeze_enabled(), freeze_position(), freeze_duration() };
    return p[index];
  }

private:
  struct
  {
    vessl::binary_p frozen;
  } params_;
  vessl::math::easing::smoother<analog_t> fader_;
  vessl::processors::delay<T> delay_;
  vessl::processors::freeze<T> freeze_;
};