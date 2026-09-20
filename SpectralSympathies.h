/**
 * Copyright 2026 Damien Quartz
 */

#pragma once 

#include "SpectralGenerator.h"
#include "BlurKernel.h"
#include "vessl/vessl.h"

template<vessl::size_t SpectrumSize, vessl::size_t Overlap = 1>
class SpectralSympathies : public vessl::unit_generator<float>, vessl::plist<6>
{
public:
  using SpectralGen = SpectralGenerator<float, SpectrumSize, Overlap>;
  using band_t  = typename SpectralGen::frequency_band;
  using SampleArray = vessl::array<float>;
  using Parameter = vessl::parameter;
  using size_t = vessl::size_t;
  using phase_t = vessl::phase_t;
  using complex_t = vessl::transform::complex<float>;
  using RippleLfo = vessl::sample::waves::unipolar::triangle<float>;
  using SmearLfo = vessl::sample::waves::bipolar::sine<float>;
  using SmearFilter = vessl::filtering::biquad<2>::high_pass<float>;
  using MeltFilter  = vessl::filtering::biquad<1>::low_pass<float>;

  static constexpr size_t overlap_size = (SpectrumSize/(Overlap*2));
  static constexpr size_t overlap_size_half = (overlap_size/2);
  static constexpr size_t smear_bands_max = 8; // SpectrumSize/128;
  
  SpectralSympathies(SpectralGen* spec_gen, float sample_rate, float* scratch_data)
    : sample_rate_(sample_rate)
    , ripple_lfo_phase_(0)
    , ripple_lfo_step_(vessl::cast<phase_t>(1.f/overlap_size))
    , smear_lfo_phase_(0)
    , smear_lfo_step_(vessl::cast<phase_t>(1.f/overlap_size))
    , filter_scratch_(scratch_data, SpectrumSize)
    , generator_(spec_gen)
  {
    params_.volume.value = 1.0f;
    params_.damping.value = 0.9f;
    params_.motion.value = 1.0f;
  }
      
  VESSL_INLINE size_t get_band_index(float frequency)
  {
    return generator_->get_band_index(frequency);
  }
  
  VESSL_INLINE float get_band_frequency(size_t band_index)
  {
    return generator_->get_band_frequency(band_index);
  }
  
  VESSL_INLINE band_t& get_band(float freq)
  {
    const size_t idx = generator_->get_band_index(freq);
    return generator_->get_band(idx);
  }

  float get_magnitude_mean()
  {
    float accum = 0;
    for (int i = 1; i < SpectrumSize/2; ++i)
    {
      accum += generator_->get_band(i).magnitude();
    }
    return (accum / (SpectrumSize/2));
  }

  VESSL_INLINE Parameter smear() const { return params_.smear("smear", 's'); }
  VESSL_INLINE Parameter damping() const { return params_.damping("damping", 'd'); }
  VESSL_INLINE Parameter melt() const { return params_.melt("melt", 'm'); }
  VESSL_INLINE Parameter ripple() const { return params_.ripple("ripple", 'r'); }
  VESSL_INLINE Parameter motion() const { return params_.motion("motion", 'o'); }
  VESSL_INLINE Parameter volume() const { return params_.volume("volume", 'v'); }

  VESSL_INLINE void excite(size_t bidx, complex_t in, float response)
  {
    //if (bidx > 1 && bidx < SpectrumSize/2)
    {
      band_t& band = generator_->get_band(bidx);

      // this has a widening effect, but is more complicated that it needs to be,
      // and can wind up overloading some bands.
      // if (band.magnitude() < response)
      // {
      //   const float in_mag = in.normalize();
      //   complex_t band_cmplx = band.to_complex();
      //   complex_t delta = in - band_cmplx;
      //   delta.scale(in_mag*response);
      //   band_cmplx.add(delta);
      //   band.set_complex(band_cmplx);
      // }

      // simple, direct, no overloading.
      if (band.magnitude() < response)
      {
        band.set_complex(in);
      }
    }
  }

  void excite(size_t bidx, float amp, phase_t phase)
  {
    if (bidx > 1 && bidx < SpectrumSize/2)
    {
      band_t& band = generator_->get_band(bidx);
      float ba = band.magnitude();
      const float ea = amp;
      if (ba < 0.8f && ea > ba)
      {
        band_t delta(ea, phase);
        delta.subtract(band);
        delta.scale(0.2f);
        band.add(delta);
      }
    }
  }
  
  void excite(float freq, float amp, phase_t phase)
  {
    const int bidx = generator_->get_band_index(freq);
    excite(bidx, amp, phase);
  }
  
  [[nodiscard]] const parameter_list& parameters() const override { return *this; }
  
  VESSL_INLINE float generate() override
  {
    // apply decay to the spectrum between overlaps.
    if (generator_->get_overlap_count() == overlap_size_half)
    {
      fill_spectrum();
    }
    
    const float volume = vessl::math::constrain(params_.volume.value, 0.f, 1.f);
    return generator_->generate()*volume;
  }
  
  VESSL_INLINE void generate(SampleArray output)
  {
    fill_spectrum();
    generator_->generate(output);
    const float volume = vessl::math::constrain(params_.volume.value, 0.f, 1.f);
    output.scale(volume);
  }

  static SpectralSympathies* create(float sample_rate)
  {
    float* scratch_data = new float[SpectrumSize];
    SpectralGen* spectral_gen = SpectralGen::create(sample_rate, vessl::sample::windows::type::triangle);
    return new SpectralSympathies(spectral_gen, sample_rate, scratch_data);
  }

  static void destroy(SpectralSympathies* synth)
  {
    SpectralGen::destroy(synth->generator_);
    delete[] synth->filter_scratch_.data();
    delete synth;
  }

protected:
  vessl::parameter element_at(vessl::size_t index) const override
  {
    switch (index)
    {
      case 0: return damping();
      case 1: return smear();
      case 2: return melt();
      case 3: return ripple();
      case 4: return motion();
      case 5: return volume();
      default: return Parameter::none();
    }
  }

private:
  VESSL_INLINE void fill_spectrum()
  {    
    const size_t count = SpectrumSize/2;
    const float smr = vessl::math::max(params_.smear.value, 0.f);
    const float mlt = params_.melt.value;
    const float ripv = vessl::math::constrain(params_.ripple.value, 0.f, 1.f);
    const float ripf = 0.5f + ripv * 1.5f;
    const float ripd = vessl::math::interp<vessl::math::easing::expo::out>(0.f, 0.025f, ripv);
    const float mot = vessl::math::max(params_.motion.value, 0.f);
    const float dmp = vessl::math::constrain(params_.damping.value + mot*0.05f + mlt*0.05f, 0.0001f, 0.9999f);

    // "melt" spectral magnitudes downwards
    for (size_t i = 1; i < count; ++i)
    {
      band_t& band = generator_->get_band(count - i);
      filter_scratch_[i-1] = band.magnitude();
    }

    const float mhz = vessl::math::lerp(sample_rate_*0.49f, sample_rate_*0.25f, mlt);
    vessl::filtering::args mlt_args(sample_rate_, mhz, vessl::filtering::q::butterworth<float>(), vessl::gain_t(0.0f));
    // apply melt
    melt_flt_.process(filter_scratch_.data(), filter_scratch_.data(), count-1, mlt_args);

    for (size_t i = 1; i < count; ++i)
    {
      band_t& band = generator_->get_band(count - i);

      // apply ripple
      const float m = filter_scratch_[i-1];
      ripple_lfo_phase_ += static_cast<phase_t>(ripple_lfo_step_*ripf);
      const phase_t rp = vessl::cast<phase_t>(m);
      const float r = ripple_lfo_.evaluate(ripple_lfo_phase_+rp);
      const float mr = vessl::math::lerp(m, r, ripd);
      band.set_magnitude(mr);
    }

    smear_lfo_phase_ += static_cast<size_t>(smear_lfo_step_*smr);
    const float smear_val = smear_lfo_.evaluate(smear_lfo_phase_);
    const float smear_mod = vessl::math::abs(smear_val)*mot;
    const bool smear_up = smear_val > 0;

    float* scratch = filter_scratch_.data();
    for (size_t i = 1; i < count; ++i)
    {
      const size_t bidx = smear_up ? i : count - i;
      band_t& band = generator_->get_band(bidx);
      complex_t cmplx = band.to_complex();
      *scratch++ = cmplx.r;
      *scratch++ = cmplx.i;
    }

    const float hz = 20.f + sample_rate_*0.25*smear_mod;
    vessl::filtering::args smr_args(sample_rate_, hz, vessl::filtering::q::butterworth<float>(), vessl::gain_t(0.0f));
    smear_flt_.process(filter_scratch_.data(), filter_scratch_.data(), (count-1)*2, smr_args);

    scratch = filter_scratch_.data();
    for (size_t i = 1; i < count; ++i)
    {
      const size_t bidx = smear_up ? i : count - i;
      band_t& band = generator_->get_band(bidx);
      const float re = *scratch++;
      const float im = *scratch++;
      // apply damping to the reconstruction
      complex_t cmplx(re*dmp, im*dmp);
      band.set_complex(cmplx);
    }

    // smear_lfo_phase_ += smear_lfo_step;
    // float smear_mod = smear_lfo_.evaluate(smear_lfo_phase_)*(smear_bands_max/2);
    // float smear_scale = vessl::math::interp<vessl::math::easing::expo::out>(8.0f, 0.125f, dmp);
    // float smear_amt = params_.spread.value * smear_scale * (1.f / Overlap);
    // const size_t smear_width = static_cast<size_t>(smear_bands_max/2 + smear_mod) * 2;
    // if (smear_width > 0 && smear_amt > 0)
    // {
    //   for (size_t i = 2 + smear_width; i < count/2 - smear_width; i+=2)
    //   {
    //     band_t& band = generator_->get_band(i);
        
    //     // "smear" the spectrum contents by blending nearby bands
    //     const size_t li = i / smear_width;
    //     const size_t hi = i * smear_width;
    //     band_t lob = li > 0 ? generator_->get_band(li) : band_t();
    //     band_t hib = hi < count ? generator_->get_band(hi) : band_t();

    //     lob.scale(smear_amt);
    //     hib.scale(smear_amt);
    //     band.add(lob);
    //     band.add(hib);

    //     // doing this nerfs the decay effect when smear is turned up.
    //     //float mag = band.magnitude();
    //     //band.scale(mag > 0.8f ? 0.8f - smear_amt*2 : 1.0f - smear_amt*2);
    //   }
    // }
  }

  struct 
  {
    vessl::analog_p damping;
    vessl::analog_p smear;
    vessl::analog_p melt;
    vessl::analog_p ripple;
    vessl::analog_p motion;
    vessl::analog_p volume;
  } params_;
  
  float sample_rate_;
  RippleLfo ripple_lfo_;
  phase_t ripple_lfo_phase_;
  phase_t ripple_lfo_step_;
  SmearLfo smear_lfo_;
  phase_t smear_lfo_phase_;
  phase_t smear_lfo_step_;
  SmearFilter smear_flt_;
  MeltFilter melt_flt_;
  SampleArray filter_scratch_;

  SpectralGen* generator_;
  bool phase_flip_ = false;
};
