#pragma once 

#include "SpectralGenerator.h"
#include "BlurKernel.h"
#include "vessl/vessl.h"

template<vessl::size_t SpectrumSize, vessl::size_t Overlap = 1>
class SpectralSympathies : public vessl::unit_generator<float>, vessl::plist<4>
{
public:
  using SpectralGen = SpectralGenerator<float, SpectrumSize, Overlap>;
  using band_t  = typename SpectralGen::frequency_band;
  using SampleArray = vessl::array<float>;
  using Parameter = vessl::parameter;
  using size_t = vessl::size_t;
  using phase_t = vessl::phase_t;
  using complex_t = vessl::transform::complex<float>;
  using SmearLfo = vessl::sample::waves::unipolar::triangle<float>;

  static constexpr size_t overlap_size = (SpectrumSize/(Overlap*2));
  static constexpr size_t overlap_size_half = (overlap_size/2);
  static constexpr size_t smear_bands_max = SpectrumSize/128;
  static constexpr phase_t smear_lfo_step = vessl::phase_180 / Overlap / smear_bands_max;
  
  SpectralSympathies(SpectralGen* spec_gen, float sample_rate)
    : sample_rate_(sample_rate)
    , smear_lfo_phase_(0)
    , generator_(spec_gen)
  {
    params_.volume.value = 1.0f;
    params_.damping.value = 0.9f;
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

  VESSL_INLINE Parameter spread() const { return params_.spread("spread", 's'); }
  VESSL_INLINE Parameter damping() const { return params_.damping("damping", 'd'); }
  VESSL_INLINE Parameter melt() const { return params_.melt("melt", 'm'); }
  VESSL_INLINE Parameter volume() const { return params_.volume("volume", 'v'); }

  VESSL_INLINE void excite(size_t bidx, complex_t in, float response)
  {
    band_t& band = generator_->get_band(bidx);
    //float band_mag = band.magnitude();
    float in_mag = in.normalize();
    //if (band_mag < 0.25f)
    {
      complex_t band_cmplx = band.to_complex();
      complex_t delta = in - band_cmplx;
      delta.scale(in_mag*response);
      band_cmplx.add(delta);
      band.set_complex(band_cmplx);
    }
  }

  void excite(size_t bidx, float amp, phase_t phase)
  {
    //if (bidx > 1 && bidx < SpectrumSize/2)
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
    if (phase_flip_)
    {
      generator_->template generate<true>(output);
      phase_flip_ = false;
    }
    else
    {
      generator_->template generate<false>(output);
      phase_flip_ = true;
    }

    const float volume = vessl::math::constrain(params_.volume.value, 0.f, 1.f);
    output.scale(volume);
  }

  static SpectralSympathies* create(float sample_rate)
  {
    SpectralGen* spectral_gen = SpectralGen::create(sample_rate, vessl::sample::windows::type::triangle);
    return new SpectralSympathies(spectral_gen, sample_rate);
  }

  static void destroy(SpectralSympathies* synth)
  {
    SpectralGen::destroy(synth->generator_);
    delete synth;
  }

protected:
  vessl::parameter element_at(vessl::size_t index) const override
  {
    switch (index)
    {
      case 0: return damping();
      case 1: return spread();
      case 2: return melt();
      case 3: return volume();
      default: return Parameter::none();
    }
  }

private:
  VESSL_INLINE void fill_spectrum()
  {    
    const float mlt = params_.melt.value*0.75f;
    const float dmp = params_.damping.value;
    const size_t count = SpectrumSize/2;
    for (size_t i = 1; i < count; ++i)
    {
      band_t& band = generator_->get_band(i);
      
      //"melt" some of this band's energy into the band below.
      const size_t mi = i == 1 ? count - 1 : i-1;
      band_t& target = generator_->get_band(mi);
      float bmag = band.magnitude();
      target.set_magnitude(target.magnitude() + bmag*mlt);
      band.set_magnitude(bmag - bmag*mlt);

      // now apply normal decay to this band
      band.scale(dmp);
    }

    smear_lfo_phase_ += smear_lfo_step;
    float smear_mod = smear_lfo_.evaluate(smear_lfo_phase_)*(smear_bands_max/2);
    float smear_scale = vessl::math::interp<vessl::math::easing::expo::out>(8.0f, 0.125f, dmp);
    float smear_amt = params_.spread.value * smear_scale * (1.f / Overlap);
    const size_t smear_width = static_cast<size_t>(smear_bands_max/2 + smear_mod) * 2;
    if (smear_width > 0 && smear_amt > 0)
    {
      for (size_t i = 2 + smear_width; i < count/2 - smear_width; i++)
      {
        band_t& band = generator_->get_band(i);
        
        // "smear" the spectrum contents by blending nearby bands
        const size_t li = i / smear_width;
        const size_t hi = i * smear_width;
        band_t lob = li > 0 ? generator_->get_band(li) : band_t();
        band_t hib = hi < count ? generator_->get_band(hi) : band_t();

        lob.scale(smear_amt);
        hib.scale(smear_amt);
        band.add(lob);
        band.add(hib);

        // doing this nerfs the decay effect when smear is turned up.
        //float mag = band.magnitude();
        //band.scale(mag > 0.8f ? 0.8f - smear_amt*2 : 1.0f - smear_amt*2);
      }
    }
  }

  struct 
  {
    vessl::analog_p damping;
    vessl::analog_p spread;
    vessl::analog_p melt;
    vessl::analog_p volume;
  } params_;
  
  float sample_rate_;
  SmearLfo smear_lfo_;
  phase_t smear_lfo_phase_;

  SpectralGen* generator_;
  bool phase_flip_ = false;
};
