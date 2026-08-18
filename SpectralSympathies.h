#pragma once 

#include "SpectralGenerator.h"
#include "BlurKernel.h"
#include "vessl/vessl.h"

template<vessl::size_t SpectrumSize, vessl::size_t Overlap = 1, bool LinearDecay = true>
class SpectralSympathies : public vessl::unit_generator<float>, vessl::plist<4>
{
public:
  using SpectralGen = SpectralGenerator<float, SpectrumSize, Overlap>;
  using band_t  = typename SpectralGen::frequency_band;

  using Spectrum = vessl::array<band_t>;
  using SampleArray = vessl::array<float>;
  using Parameter = vessl::parameter;
  using size_t = vessl::size_t;
  using phase_t = vessl::phase_t;
  using complex_t = vessl::transform::complex<float>;
  
  // all data arrays should be at least bands_size long.
  SpectralSympathies(SpectralGen* spec_gen, float sample_rate)
    : sample_rate_(sample_rate)
    , generator_(spec_gen)
    , overlap_size_(SpectrumSize/(Overlap*2))
    , overlap_size_half_(overlap_size_/2)
  {
    params_.volume.value = 1.0f;
    params_.decay.value = vessl::duration_t::from_seconds(1.0f, sample_rate);
    set_decay(1.0f);
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
  VESSL_INLINE Parameter decay() const { return params_.decay("decay", 'd'); }\
  VESSL_INLINE Parameter melt() const { return params_.melt("melt", 'm'); }
  VESSL_INLINE Parameter volume() const { return params_.volume("volume", 'v'); }

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
    float decay_param = params_.decay.value.to_seconds(sample_rate_);
    if (vessl::math::abs(decay_seconds_ - decay_param) > 0.001f)
    {
      set_decay(decay_param);
    }
    
    // apply decay to the spectrum between overlaps.
    if (generator_->get_overlap_count() == overlap_size_half_)
    {
      fill_spectrum();
    }
    
    const float volume = vessl::math::constrain(params_.volume.value, 0.f, 1.f);
    return generator_->generate()*volume;
  }
  
  VESSL_INLINE void generate(SampleArray output)
  {
    auto w = output.make_writer();
    while (w)
    {
      w << generate();
    }
  }

  static SpectralSympathies* create(float sample_rate)
  {
    SpectralGen* spectral_gen = SpectralGen::create(sample_rate, vessl::sample::windows::type::triangle);
    return new SpectralSympathies(spectral_gen, sample_rate);
  }

  static void destroy(SpectralSympathies* synth)
  {
    SpectralGen::destroy(synth->generator_);
    delete[] synth->spectrum_.data();
    delete synth;
  }

protected:
  vessl::parameter element_at(vessl::size_t index) const override
  {
    switch (index)
    {
      case 0: return decay();
      case 1: return spread();
      case 2: return melt();
      case 3: return volume();
      default: return Parameter::none();
    }
  }

private:
  void set_decay(const float in_seconds)
  {
    // having a shorter decay than the overlap size doesn't make sense
    // and we also want to avoid divide-by-zero.
    decay_seconds_ = vessl::math::max(overlap_size_ / sample_rate_, in_seconds);
    if constexpr (LinearDecay)
    {
      // amplitude needs to decrease by 1 / (decaySeconds * sampleRate()) every sample.
      // eg decaySeconds == 1 -> 1 / sampleRate()
      //    decaySeconds == 0.5 -> 1 / (0.5 * sampleRate), which is twice as fast, equivalent to 2 / sampleRate()
      // since we generate a new buffer every overlapSize samples, we multiply that rate by overlapSize, giving:
      decay_dec_ = overlap_size_ / (decay_seconds_ * sample_rate_);
    }
    else // exponential decay
    {
      float block_rate = sample_rate_ / overlap_size_;
      float length_in_blocks = decay_seconds_ * block_rate;
      decay_dec_ = 1.0 + vessl::math::log(0.0001f) / (length_in_blocks + 20);
    }
  }

  VESSL_INLINE void fill_spectrum()
  {    
    const float mlt = params_.melt.value;
    const size_t count = SpectrumSize/2;
    for (size_t i = 1; i < count; ++i)
    {
      band_t& band = generator_->get_band(i);
      
      // "melt" some of this band's energy into the band below,
      // wrapping around to the top of the spectrum if we are at the bottom.
      const size_t mi = i == 1 ? count - 1 : i-1;
      band_t& target = generator_->get_band(mi);
      band_t  delt = band;
      delt.scale(mlt);
      band.scale(1.0f - mlt);
      target.add(delt);

      // now apply normal decay to this band
      band.scale(decay_dec_);
    }

    // @todo can overwhelm the sound pretty easily
    const size_t smr = static_cast<size_t>(params_.spread.value*32) * 2;
    if (smr > 0)
    {
      for (size_t i = 2 + smr; i < count/2 - smr; i++)
      {
        band_t& band = generator_->get_band(i);
        
        // "smear" the spectrum contents by blending nearby bands
        const size_t li = i / smr;
        const size_t hi = i * smr;
        band_t lob = li > 0 ? generator_->get_band(li) : band_t();
        band_t hib = hi < count ? generator_->get_band(hi) : band_t();

        lob.scale(0.2f);
        hib.scale(0.2f);
        band.add(lob);
        band.add(hib);
        float mag = band.magnitude();
        band.scale(mag > 0.8f ? 0.4f : 0.6f);
      }
    }
  }

  // VESSL_INLINE void process_band(int idx, int spec_size)
  // {
  //   Band& b = bands_[idx];
  //   if (LinearDecay)
  //   {
  //     b.decay = b.decay > decay_dec_ ? b.decay - decay_dec_ : 0;
  //   }
  //   else
  //   {
  //     //b.decay *= decayDec;
  //     b.amplitude *= decay_dec_;
  //   }
  //   
  //   //if (b.decay > 0)
  //   {
  //     float bright = params_.brightness.value;
  //     //float a = b.decay*b.amplitude;
  //     float a = b.amplitude;
  //     spec_bright_[idx] += a;
  //     constexpr int iters = SpectralBandPartials;
  //     for (int i = 0; i < iters && b.partials[i] < spec_size; ++i)
  //     {
  //       int p = 2 + i;
  //       a *= bright;
  //       int pidx = b.partials[i];
  //       spec_bright_[pidx] += a / p;
  //     }
  //   }
  // }

  struct 
  {
    vessl::duration_p decay;
    vessl::analog_p   spread;
    vessl::analog_p   melt;
    vessl::analog_p   volume;
  } params_;
  
  float sample_rate_;
  // cache this so we only recalculate decay_dec_ when necessary.
  float decay_seconds_;
  float decay_dec_;

  // cached so we only recalc the kernel when needed
  float smear_amount_;
  
  Spectrum     spectrum_;
  SpectralGen* generator_;
  
  size_t overlap_size_;
  size_t overlap_size_half_;
  size_t overlap_size_mask_;
};
