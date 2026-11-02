#pragma once

#include "vessl/vessl.h"
#include "vessl/vessl_units.h"

template<typename T>
class ClockedDelay : public vessl::unit_processor<T>
                   , public vessl::time::clockable
                   , protected vessl::plist<2>
{
  using analog_t = vessl::analog_t;
  using Array = vessl::array<T>;
  using Delay = vessl::processors::delay<T>;
  using Parameter = vessl::parameter;

public:
  ClockedDelay(Array delay_buffer, analog_t sample_rate, period_t clock_period_min, period_t clock_period_max, analog_t bpm = 60, analog_t feedback_amount = 0, vessl::unit_processor<T>* feedback_proc = nullptr)
  : clockable(sample_rate, clock_period_min, clock_period_max, bpm)
  , delay_(delay_buffer, sample_rate_, tempo_.to_seconds(sample_rate), feedback_amount, feedback_proc)
  {
    params_.ratio.value = 1.f;
  }
  
  [[nodiscard]] const vessl::parameter_list& parameters() const override { return *this; }

  // set delay time as ratio of the bpm
  [[nodiscard]] VESSL_INLINE Parameter ratio() const { return params_.ratio("ratio", 'r'); }
  [[nodiscard]] VESSL_INLINE Parameter feedback() const { return delay_.feedback(); }

  T process(const T& in) override
  {
    tick();
    delay_.time() = tempo_.samples * params_.ratio.value;
    return delay_.process(in);
  }

  VESSL_INLINE void process(vessl::source<T>& source, vessl::sink<T>& sink) override { vessl::processor<T>::process(source, sink); }
    
  template<vessl::time::mode TimeMode>
  void process(Array input, Array output)
  {
    tick(input.size());
    delay_.time() = tempo_.samples * params_.ratio.value;
    delay_.template process<TimeMode>(input, output);
  }

protected:
  [[nodiscard]] Parameter element_at(vessl::size_t index) const override
  {
    switch(index)
    {
      case 0: return feedback();
    }

    return Parameter::none();
  }

private:
  struct
  {
    vessl::analog_p ratio;
  } params_;

  Delay delay_;
};