#pragma once

#include "vessl/vessl.h"
#include "Program.h"


// Evaluator operates on sample frames
// because the underlying Program requires this.
template<typename T, vessl::size_t N>
class Evaluator : public vessl::unit_processor<vessl::sample::frame<T,N>>, public vessl::plist<11>
{
public:
  using size_t   = vessl::size_t;
  using analog_t = vessl::analog_t;
  using input_t  = vessl::sample::frame<T,N>;
  using output_t = vessl::sample::frame<T,N>;
  using Parameter = vessl::parameter;
  using Smoother = vessl::math::easing::smoother<float>;
  using Exec = Program::Executable<512, 4096>;

  Evaluator(const analog_t sample_rate)
  : sample_rate_(sample_rate)
  , tick_phase_(0)
  {
    params_.tempo.value = 120;
    params_.bit_depth.value = 16;
    params_.rate.value = 1;
  }

  ~Evaluator()
  {
  }

  [[nodiscard]] Parameter tick() const { return params_.tick("Tick", 't'); }
  [[nodiscard]] Parameter tempo() const { return params_.tempo("Tempo", 'q'); }
  [[nodiscard]] Parameter bits() const { return params_.bit_depth("Bits", 'w'); }
  [[nodiscard]] Parameter rate() const { return params_.rate("Rate", 'r'); }
  [[nodiscard]] Parameter v(size_t idx) const { return idx < 8 ? params_.v[idx]("V", '0' + idx) : Parameter::none(); }

  [[nodiscard]] const vessl::parameter_list& parameters() const override { return *this; }

  // return true if compilation succeeds.
  // on failure can retrive error code and error position from getters.
  bool compile(const char* program_text)
  {
    std::vector<Program::Op> ops;
    compile_result_ = Program::CE_NONE;
    compile_error_position_ = -1;
    Program::Compile<512, 4096>(program_text, ops, compile_result_, compile_error_position_);
    bool success = compile_result_ == Program::CE_NONE;
    if (success)
    {
      exec_.Load(ops, vessl::math::random::u32());
      params_.tick.value = 0;
      tick_phase_ = 1.f;
    }
    return success;
  } 

  [[nodiscard]] Program::CompileError last_compile_error() const { return compile_result_; }
  [[nodiscard]] int last_compile_error_position() const { return compile_error_position_; }
  [[nodiscard]] Program::RuntimeError last_runtime_error() const { return run_result_; }
  [[nodiscard]] Program::Value get(Program::Char var) const { return exec_.Get(var); }
  [[nodiscard]] size_t instruction_count() const { return exec_.GetInstructionCount(); }

  static Program::Compiler* make_compiler(const char* program_text)
  {
    Program::Compiler* c = Program::Compile<512, 4096>(program_text);
    return c;
  }

  bool load(Program::Compiler* result)
  {
    compile_result_ = result->ErrorCode();
    compile_error_position_ = result->ErrorLocation();
    bool success = compile_result_ == Program::CE_NONE;
    if (success)
    {
      exec_.Load(result->Ops(), vessl::math::random::u32());
      params_.tick.value = 0;
      tick_phase_ = 1.f;
    }
    return success;
  }

  VESSL_INLINE output_t process(const input_t& in) override
  {
    output_t out(0);
    { 
      Program::Value range = static_cast<Program::Value>(vessl::math::pow(2.f, params_.bit_depth.value));
      float inv_range = 1.0f / range;

      tick_phase_ += params_.rate.value;
      if (tick_phase_ >= 1.f)
      {
        Program::Value tick_adv = static_cast<Program::Value>(tick_phase_);
        Program::Value tick_now = params_.tick.value + tick_adv;
        tick_phase_ -= static_cast<analog_t>(tick_adv);

        for(size_t i = 0; i < N; ++i)
        {
          float norm = vessl::cast<float>(in.samples[i])*0.5f + 0.5f;
          args_[i] = static_cast<Program::Value>(norm * (range-1));
        }
        const float sr = sample_rate_;
        const float mdenom = 1000.f / sr;
        const float qdenom = 128.f / (sr / (params_.tempo.value / 60.f));
        exec_['t'] = params_.tick.value;
        exec_['m'] = (Program::Value)vessl::math::round(params_.tick.value * mdenom);
        exec_['q'] = (Program::Value)vessl::math::round(params_.tick.value * qdenom);
        exec_['w'] = range;
        exec_['~'] = (Program::Value)sr;
        
        for(size_t vi = 0; vi < 8; ++vi)
        {
          exec_.SetVC(vi, static_cast<Program::Value>(255.f * params_.v[vi].value));
        }

        run_result_ = exec_.Run(args_, N);
        params_.tick.value = tick_now;
      }

      for(size_t i = 0; i < N; ++i)
      {
        float norm = (args_[i] % range) * inv_range; 
        out.samples[i] = vessl::cast<T>(norm*2.f - 1.f);
      }
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
      default: return v(index - 3);
    }

    return Parameter::none();
  }

private:
  struct 
  {
    vessl::param<Program::Value> tick;
    vessl::analog_p tempo;
    vessl::analog_p bit_depth;
    vessl::analog_p rate;
    vessl::analog_p v[8];
  } params_;

  analog_t sample_rate_;
  analog_t tick_phase_;
  Exec     exec_;
  Program::Value args_[N];
  Program::CompileError compile_result_;
  Program::RuntimeError run_result_;
  int compile_error_position_;
};