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
    Synthesizes sound by using overlap-add IFFT synthesis of a spectrum
    that is "excited" by the spectrum derived from the analyzed input.
*/

#pragma once

#include "vessl/vessl.h"
#include "SpectralSympathies.h"

// @todo - so a thing that sounds pretty cool is reducing the SpectrumSize down to like 512.
// It creates more of a talkbox kind of effect.
// What I want to try is:
// Density is a blend between a set of small forward FFTs.
// The generator continues to run at 4096, but ideally with higher overlap.
// Or dynamic overlap based on Density?
// Equivalent might be to group the forward FFT bands into average bands, which are our strings.
// So we can run the forward FFT at the same size as the generator,
// but when exciting it we are using average band information rather than a single band's information.
// The string mapping would then blend from matching center frequency to center frequency
// to a linear mapping from string index to generator band index.
//
// Another thot: spread could be around the strings from the source FFT, so that when exciting
// adjacent bands we are doing so with real data.
// The effect the SpectralSympathies generator is doing with spread is more like a "smear".

template<typename T, size_t SpectrumSize>
class Condolences : public vessl::unit_processor<T>, public vessl::plist<6>
{
public:
  using size_t = vessl::size_t;
  using analog_t = vessl::analog_t;
  using sample_t = T;
  using complex_t = vessl::transform::complex<sample_t>;
  using Smoother = vessl::math::easing::smoother<vessl::analog_t>;
  using Parameter = vessl::parameter;
  using SpectralGen = SpectralSympathies<SpectrumSize, false>;
  using SampleArray = vessl::array<sample_t>;
  using ComplexArray = vessl::array<complex_t>;
  using Window = vessl::sample::windows::type;
  using FFT = vessl::transform::fft<sample_t>;
  using Frequency = vessl::frequency<analog_t>;

  static constexpr size_t spread_width = 4;
  static constexpr float spread_pct = 0.1f;
  
  Condolences(
    const analog_t sample_rate,
    const size_t block_size, 
    sample_t* input_buffer_data, 
    sample_t* input_window_data,
    sample_t* input_analyze_data, 
    complex_t* input_spectrum_data, 
    complex_t* feedback_spectrum_data,
    SpectralGen* spectral_generator
  )
  : sample_rate_(sample_rate)
  , density_min_(16)
  , density_max_(static_cast<float>(SpectrumSize)/16.f)
  , band_first_idx_(1.f + spread_width)
  , band_last_idx_(static_cast<float>(SpectrumSize/2) - spread_width - 1)
  , decay_min_(static_cast<float>(SpectrumSize) * 0.5f / sample_rate)
  , input_buffer_write_(0)
  , input_buffer_(input_buffer_data, SpectrumSize)
  , input_window_(input_window_data, SpectrumSize)
  , input_analyze_(input_analyze_data, SpectrumSize)
  , input_spectrum_(input_spectrum_data, SpectrumSize / 2)
  , feedback_spectrum_(feedback_spectrum_data, SpectrumSize / 2)
  , input_transform_(SpectrumSize)
  , spectral_gen_(spectral_generator)
  {
    vessl::sample::windows::render(Window::hann, input_window_);
  }

  [[nodiscard]] Parameter density() const { return params_.density("density", 'd'); }
  [[nodiscard]] Parameter spacing() const { return params_.spacing("spacing", 's'); }
  [[nodiscard]] Parameter spread() const { return params_.spread("spread", 'r'); }
  [[nodiscard]] Parameter melt() const { return params_.melt("melt", 'm'); }
  [[nodiscard]] Parameter decay() const { return params_.decay("decay", 'c'); }
  [[nodiscard]] Parameter feedback() const { return params_.feedback("feedback", 'f'); }
  
  [[nodiscard]] const parameter_list& parameters() const override { return *this; }
  
  VESSL_INLINE sample_t process(const sample_t& in) override
  {
    VASSERT(false, "Condolences only supports block processing arrays.");
    return in;
  }
  
  VESSL_INLINE void process(vessl::array<T> in, vessl::array<T> out) override
  {
    const size_t block_size = in.size();
    
    density_ = vessl::math::lerp(density_min_,  density_max_, params_.density.value);
    spacing_ = params_.spacing.value;
    spread_ = vessl::math::interp<vessl::math::easing::quad::out>(0.f, 1.f, params_.spread.value);
    spread_max_ = vessl::math::lerp(SpectrumSize/4.f, SpectrumSize/64.f, params_.density.value);
    decay_ = vessl::math::max(decay_min_, params_.decay.value);
    melt_ = params_.melt.value;
    
    // reduce volume based on combination of decay, spread, and brightness parameters
    volume_ = vessl::math::interp<vessl::math::easing::expo::out>(1.0f, 0.5f, 
        0.2f*params_.decay.value
      + 0.2f*params_.spread.value);
    
    //spectral_gen_->spread() = spread_.value;
    //spectral_gen_->set_spread_bands_max(spread_max_.value);
    spectral_gen_->decay() = vessl::duration_t::from_seconds(decay_.value, sample_rate_);
    spectral_gen_->melt() = melt_.value;
    spectral_gen_->volume() = volume_.value;
    
    const size_t string_count = vessl::math::max(get_string_count(), 1ull);
    constexpr analog_t mag_norm = 256.f / static_cast<float>(SpectrumSize);
    for (size_t i = 0; i < block_size; ++i)
    {
      input_buffer_[input_buffer_write_++] = in[i];
      if (input_buffer_write_ == SpectrumSize)
      {
        // window the input and output to an analysis buffer
        // because running the fft messes up the input samples.
        input_window_.multiply(input_buffer_, input_analyze_);
        input_transform_.forward(input_analyze_, input_spectrum_);
        
        // transfer spectrum data from input analysis to spectral_gen
        // by sampling only those frequencies represented by our strings.
        // i.e. comb filter it.
        const size_t iss = input_spectrum_.size();
        for (size_t si = 0; si < string_count; ++si)
        {
          const float freq = frequency_of_string(si);
          const size_t bi = spectral_gen_->get_band_index(freq);
          if (bi > 0 && bi < iss)
          {
            float mag = mag_norm;
            complex_t input = input_spectrum_[bi];
            float in_mag = input.magnitude() * mag;
            vessl::phase_t in_phase = input.phase();
            spectral_gen_->excite(bi, in_mag, in_phase);
            for(size_t si = 1; si < spread_width + 1; ++si)
            {
              mag *= spread_.value;
              size_t hi = bi+si;
              size_t lo = bi-si;
              
              input = input_spectrum_[lo];
              spectral_gen_->excite(lo, input.magnitude()*mag, input.phase());

              input = input_spectrum_[hi];
              spectral_gen_->excite(hi, input.magnitude()*mag, input.phase());
            }
          }
        }
        
        // copy the back half of the array to the front half
        // continue recording input from the middle of the array.
        // doing this means we can update the spectral data for sound generation every overlap.
        input_buffer_write_ = SpectrumSize / 2;
        SampleArray input_buffer_back(input_buffer_.data() + input_buffer_write_, input_buffer_write_);
        input_buffer_back.copy_to(input_buffer_);
      }
      out[i] = spectral_gen_->generate();
    }
  }

  VESSL_INLINE analog_t get_decay_min() const { return decay_min_; };
  
  VESSL_INLINE typename SpectralGen::band_t get_band(analog_t freq_in_hz) const
  {
    return spectral_gen_->get_band(freq_in_hz);
  }
  
  // get the current string count based on the density setting
  VESSL_INLINE size_t get_string_count() const
  {
    return static_cast<size_t>(density_.value + 0.5f);
  }

  VESSL_INLINE analog_t frequency_of_string(const int string_num) const
  {
    const float t = static_cast<float>(string_num) / get_string_count();
    // convert first and last bands to midi notes and then do a linear interp, converting back to Hz at the end.
    const Frequency low_freq = Frequency::of_hertz(spectral_gen_->get_band_frequency(band_first_idx_), sample_rate_);
    const Frequency hi_freq = Frequency::of_hertz(spectral_gen_->get_band_frequency(band_last_idx_), sample_rate_);
    const float lin_freq = vessl::math::lerp(low_freq.norm, hi_freq.norm, t);
    // @todo given normalized frequency, we probably don't need to convert to midi note to get log frequency
    const float midi_note = vessl::math::lerp(low_freq.as_midi_note(sample_rate_), hi_freq.as_midi_note(sample_rate_), t);
    const float log_freq = Frequency::of_midi_note(midi_note, sample_rate_).norm;
    // we lerp from logFreq up to linFreq because log spacing clusters frequencies
    // towards the bottom of the range, which means that when holding down the mouse on a string
    // and lowering this param, you'll hear the pitch drop, which makes more sense than vice-versa.
    return vessl::math::lerp(log_freq, lin_freq, spacing_.value) * sample_rate_;
  }
  
  static Condolences* create(vessl::analog_t sample_rate, vessl::size_t block_size)
  {
    sample_t* input_buffer_data = new sample_t[SpectrumSize];
    sample_t* input_window_data = new sample_t[SpectrumSize];
    sample_t* input_analyze_data = new sample_t[SpectrumSize];
    complex_t* input_spectrum_data = new complex_t[SpectrumSize/2];
    complex_t* feedback_spectrum_data = new complex_t[SpectrumSize/2];
    SpectralGen* spectral_generator = SpectralGen::create(sample_rate);
    return new Condolences(sample_rate, block_size,
      input_buffer_data, 
      input_window_data, 
      input_analyze_data,
      input_spectrum_data,
      feedback_spectrum_data,
      spectral_generator
      );
  }
  
  static void destroy(const Condolences* condolences)
  {
    if (condolences)
    {
      SpectralGen::destroy(condolences->spectral_gen_);
      delete[] condolences->feedback_spectrum_.data();
      delete[] condolences->input_spectrum_.data();
      delete[] condolences->input_analyze_.data();
      delete[] condolences->input_window_.data();
      delete[] condolences->input_buffer_.data();
    }
  }
  
protected:
  [[nodiscard]] vessl::parameter element_at(vessl::size_t index) const override
  {
    switch (index)
    {
      case 0: return density();
      case 1: return spacing();
      case 2: return spread();
      case 3: return melt();
      case 4: return decay();
      case 5: return feedback();
      default: return Parameter::none();
    }
  }

private:
  struct
  {
    vessl::analog_p density;
    vessl::analog_p spacing;
    vessl::analog_p spread;
    vessl::analog_p melt;
    vessl::analog_p decay;
    vessl::analog_p feedback;
  } params_;
  
  Smoother density_;
  Smoother spacing_;
  Smoother spread_;
  Smoother spread_max_;
  Smoother melt_;
  Smoother decay_;
  Smoother feedback_;
  Smoother volume_;
  
  analog_t sample_rate_;
  analog_t density_min_;
  analog_t density_max_;
  analog_t band_first_idx_;
  analog_t band_last_idx_;
  analog_t decay_min_;
  
  size_t input_buffer_write_;
  SampleArray  input_buffer_;
  SampleArray  input_window_;
  SampleArray  input_analyze_;
  ComplexArray input_spectrum_;
  ComplexArray feedback_spectrum_;
  
  FFT input_transform_;

  SpectralGen* spectral_gen_;
};