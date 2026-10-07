#pragma once

#include "vessl/vessl.h"
#include "Program.h"


// Evaluator operates on sample frames
// because the underlying Program requires this.
template<typename T, vessl::size_t N>
class Evaluator : public vessl::unit_processor<vessl::sample::frame<T,N>>, public vessl::plist<3>
{
public:
  using size_t   = vessl::size_t;
  using analog_t = vessl::analog_t;
  using input_t  = vessl::sample::frame<T,N>;
  using output_t = vessl::sample::frame<T,N>;
  using Parameter = vessl::parameter;

  Evaluator(const analog_t sample_rate)
  : sample_rate_(sample_rate)
  , tick_phase_(0)
  , program_(nullptr)
  {
    params_.tempo.value = 120;
    params_.bit_depth.value = 16;
    params_.rate.value = 1;
  }

  ~Evaluator()
  {
    delete program_;
  }

  [[nodiscard]] Parameter tick() const { return params_.tick("Tick", 't'); }
  [[nodiscard]] Parameter tempo() const { return params_.tempo("Tempo", 'q'); }
  [[nodiscard]] Parameter bits() const { return params_.bit_depth("Bits", 'w'); }
  [[nodiscard]] Parameter rate() const { return params_.rate("Rate", 'r'); }

  [[nodiscard]] const vessl::parameter_list& parameters() const override { return *this; }

  // return true if compilation succeeds.
  // on failure can retrive error code and error position from getters.
  bool compile(const char* program_text, size_t memory_size)
  {
    delete program_;
    compile_result_ = Program::CE_NONE;
    compile_error_position_ = -1;
    program_ = Program::Compile(program_text, memory_size, compile_result_, compile_error_position_);
    bool success = compile_result_ == Program::CE_NONE;
    if (success)
    {
      params_.tick.value = 0;
    }
    return success;
  }

  [[nodiscard]] Program::CompileError last_compile_error() const { return compile_result_; }
  [[nodiscard]] int last_compile_error_position() const { return compile_error_position_; }

  output_t process(const input_t& in) override
  {
    output_t out(0);
    if (program_)
    {
      const Program::Value range = (Program::Value)1 << params_.bit_depth.value;
      const double mdenom = sample_rate_ / 1000.0;
      const double qdenom = (sample_rate_ / (params_.tempo.value / 60.0)) / 128.0;
      program_->Set('t', params_.tick.value);
      program_->Set('m', (Program::Value)vessl::math::round(params_.tick.value / mdenom));
      program_->Set('q', (Program::Value)vessl::math::round(params_.tick.value / qdenom));
      program_->Set('w', range);
      program_->Set('~', (Program::Value)sample_rate_);
      for(size_t i = 0; i < N; ++i)
      {
        // @todo this assumes input is floating-point, need better conversion
        program_args_[i] = (Program::Value)((in.samples[i] + 1) * (range / 2));
      }
      run_result_ = program_->Run(program_args_, 2);
      for(size_t i = 0; i < N; ++i)
      {
        out.samples[i] = (-1.0 + 2.0*((double)(program_args_[i] % range) / (range - 1)));
      }
      tick_phase_ += params_.rate.value;
      Program::Value tick_adv = static_cast<Program::Value>(tick_phase_);
      params_.tick.value += tick_adv;
      tick_phase_ -= static_cast<analog_t>(tick_adv);
    }
    return out;
  }

protected:
  [[nodiscard]] Parameter element_at(size_t index) const override
  {
    switch(index)
    {
      case 0: return tick();
      case 1: return tempo();
      case 2: return bits();
      case 3: return rate();
    }

    return Parameter::none();
  }

private:
  struct 
  {
    vessl::param<Program::Value> tick;
    vessl::analog_p tempo;
    vessl::param<uint8_t> bit_depth;
    vessl::analog_p rate;
  } params_;

  analog_t sample_rate_;
  analog_t tick_phase_;
  Program* program_;
  Program::Value program_args_[N];
  Program::CompileError compile_result_;
  Program::RuntimeError run_result_;
  int compile_error_position_;
};