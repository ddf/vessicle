#pragma once

#include "vessl/vessl.h"

// boring-ass name, but can't think of anything clever right now
template<typename T>
class Counter : public vessl::unit_generator<T>
              , public vessl::time::clockable
              , public vessl::plist<3>
{
public:
  using Parameter = vessl::parameter;
  using size_t = vessl::size_t;
  using analog_t = vessl::analog_t;

  Counter() 
  : clockable(1, 1, clockable::period_max)
  , value_(0) 
  {
    params_.inc.value = 1;
    params_.max.value = 16;
  }

  T generate() override
  {
    tick();
    if (params_.carry.value > 0)
    {
      params_.carry.value -= 1;
    }
    return value_;
  };

  [[nodiscard]] Parameter increment() const { return params_.inc("Increment", 'i'); }
  [[nodiscard]] Parameter maximum() const { return params_.max("Maximum", 'm'); }
  [[nodiscard]] Parameter carry() const { return params_.carry("Carry", 'c'); }
  [[nodiscard]] T value() const { return value_; }

  [[nodiscard]] const parameter_list& parameters() const { return *this; }

protected:
  void tock(period_t sample_delay)
  {
    value_ += params_.inc.value;
    if (value_ > params_.max.value)
    {
      params_.carry.value = value_ - params_.max.value;
      value_ = value_ % params_.max.value;
    }
  }

  [[nodiscard]] Parameter element_at(size_t index) const override
  {
    switch(index)
    {
      case 0: return increment();
      case 1: return maximum();
      case 2: return carry();
    }

    return Parameter::none();
  }

private:
  struct 
  {
    vessl::param<T> inc;
    vessl::param<T> max;
    vessl::param<T> carry;
  } params_;

  T value_;
};