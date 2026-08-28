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
    that is "excited" by the spectrum derived from the analyzed input
    and then modified in various ways.
*/

#pragma once

#include "vessl/vessl.h"
#include "SpectralSympathies.h"

/** @todo
  - revisit volume adjustment formula
*/

template<typename T, uint16_t SpectrumSize, uint8_t Overlap>
class Condolences : public vessl::unit_processor<T>, public vessl::plist<7>
{
public:
  using sample_t     = T;
  using size_t       = vessl::size_t;
  using analog_t     = vessl::analog_t;
  using complex_t    = vessl::transform::complex<sample_t>;
  using Smoother     = vessl::math::easing::smoother<vessl::analog_t>;
  using Parameter    = vessl::parameter;
  using Sympathies   = SpectralSympathies<SpectrumSize, Overlap>;
  using SampleArray  = vessl::array<sample_t>;
  using ComplexArray = vessl::array<complex_t>;
  using StringArray  = vessl::array<uint16_t>;
  using Window       = vessl::sample::windows::type;
  using FFT          = vessl::transform::fft<sample_t>;
  using Frequency    = vessl::frequency<analog_t>;

  static constexpr size_t AnalysisSize = Sympathies::overlap_size;
  static constexpr size_t StringCountMax = AnalysisSize/2;
  static constexpr size_t SpreadWidth = 2;
  // for clamping the param
  static constexpr float DensityMin = 4;
  static constexpr float DensityMax = static_cast<float>(AnalysisSize/2)/(SpreadWidth*2);
  
  Condolences(
    const analog_t sample_rate,  
    sample_t* input_window_data,
    sample_t* input_buffer_data,
    complex_t* input_spectrum_data,
    Sympathies* spectral_generator
  )
  : sample_rate_(sample_rate)
  , band_spacing_(sample_rate/AnalysisSize)
  , band_first_idx_(1.f + SpreadWidth)
  , band_last_idx_(static_cast<float>(AnalysisSize/2) - SpreadWidth - 1)
  , decay_min_(static_cast<float>(Sympathies::overlap_size) / sample_rate)
  , input_window_(input_window_data, AnalysisSize)
  , input_buffer_(input_buffer_data, AnalysisSize)
  , input_spectrum_(input_spectrum_data, AnalysisSize/2)
  , input_fft_(AnalysisSize)
  , spectral_gen_(spectral_generator)
  , input_buffer_write_idx_(0)
  {
    params_.response.value = 0.1f;
    input_buffer_.fill(0);
  }

  [[nodiscard]] Parameter density() const { return params_.density("density", 'd'); }
  [[nodiscard]] Parameter spacing() const { return params_.spacing("spacing", 's'); }
  [[nodiscard]] Parameter spread() const { return params_.spread("spread", 'r'); }
  [[nodiscard]] Parameter smear() const { return params_.smear("smear", 'e'); }
  [[nodiscard]] Parameter melt() const { return params_.melt("melt", 'm'); }
  [[nodiscard]] Parameter sensitivity() const { return params_.response("sensitivity", 't'); }
  [[nodiscard]] Parameter decay() const { return params_.decay("decay", 'c'); }
  
  [[nodiscard]] const parameter_list& parameters() const override { return *this; }

  [[nodiscard]] float get_input_band_magnitude(float band_freq) const
  {
    const size_t bidx = spectral_gen_->get_band_index(band_freq);
    return input_spectrum_[bidx].magnitude();
  }
  
  VESSL_INLINE sample_t process(const sample_t& in) override
  {
    VASSERT(false, "Condolences only supports block processing arrays.");
    return in;
  }
  
  VESSL_INLINE void process(vessl::array<T> in, vessl::array<T> out) override
  {
    const size_t block_size = in.size();
    
    smear_    = params_.smear.value;
    spread_   = vessl::math::interp<vessl::math::easing::quad::out>(0.f, 1.f, params_.spread.value);
    response_ = params_.response.value;
    melt_     = params_.melt.value;
    damping_  = get_damping(params_.decay.value, 0.f);

    density_ = vessl::math::constrain(params_.density.value, DensityMin, DensityMax);
    spacing_ = params_.spacing.value;
    
    // reduce volume based on combination of decay, spread, and brightness parameters
    // volume_ = vessl::math::interp<vessl::math::easing::expo::out>(1.0f, 0.5f, 
    //     0.2f*params_.decay.value
    //   + 0.2f*params_.spread.value);
    volume_ = 0.25f;
    
    spectral_gen_->spread() = smear_.value;
    spectral_gen_->damping() = damping_.value;
    spectral_gen_->melt() = melt_.value;
    spectral_gen_->volume() = volume_.value;

    const size_t string_count = vessl::math::max(static_cast<size_t>(density_.value), 1ull);
    for(size_t i = 0; i < block_size; ++i)
    {
      out[i] = input_buffer_[input_buffer_write_idx_];
      input_buffer_[input_buffer_write_idx_] = in[i]*input_window_[input_buffer_write_idx_];
      ++input_buffer_write_idx_;

      if (input_buffer_write_idx_ == AnalysisSize)
      {
        input_fft_.forward(input_buffer_, input_spectrum_);

        // transfer spectrum data from input analysis to spectral_gen
        // by sampling only those frequencies represented by our strings.
        // i.e. comb filter it.
        size_t si = 0;
        uint16_t pbi = 0;
        while(si < string_count)
        {
          const float st = static_cast<float>(si)/(string_count-1);
          const size_t abi = static_cast<size_t>(
            vessl::math::interp<vessl::math::easing::expo::in>(band_first_idx_, band_last_idx_, st)
          );

          const size_t fbi = abi > pbi ? abi : pbi+1;
          float response = response_.value;
          // main string
          {
            const size_t tbi = spectral_gen_->get_band_index(fbi*band_spacing_); 
            spectral_gen_->excite(tbi, input_spectrum_[fbi], response);
          }

          // spread strings
          {
            response *= spread_.value;
            const size_t tbi0 = f2t(fbi-1);
            const size_t tbi1 = f2t(fbi+1);
            spectral_gen_->excite(tbi0, input_spectrum_[fbi-1], response);
            spectral_gen_->excite(tbi1, input_spectrum_[fbi+1], response);
          }
          {
            response *= spread_.value;
            const size_t tbi0 = f2t(fbi-2);
            const size_t tbi1 = f2t(fbi+2);
            spectral_gen_->excite(tbi0, input_spectrum_[fbi-2], response);
            spectral_gen_->excite(tbi1, input_spectrum_[fbi+2], response);
          }

          ++si;
          pbi = fbi;
        }

        spectral_gen_->generate(input_buffer_);
        input_buffer_write_idx_ = 0;
      }
    }
  }
  
  VESSL_INLINE typename Sympathies::band_t get_band(analog_t freq_in_hz) const
  {
    return spectral_gen_->get_band(freq_in_hz);
  }
  
  static Condolences* create(vessl::analog_t sample_rate, Window input_window_type = Window::hann)
  {
    // allocate the generator first because it needs the largest contiguous block of memory
    Sympathies* spectral_generator  = Sympathies::create(sample_rate);
    sample_t*   input_window_data   = new sample_t[AnalysisSize];
    sample_t*   input_buffer_data   = new sample_t[AnalysisSize];
    complex_t*  input_spectrum_data = new complex_t[AnalysisSize/2];
    
    vessl::sample::windows::render(input_window_type, input_window_data, AnalysisSize);
    
    return new Condolences(sample_rate,
      input_window_data,
      input_buffer_data, 
      input_spectrum_data,
      spectral_generator
      );
  }
  
  static void destroy(const Condolences* condolences)
  {
    if (condolences)
    {
      Sympathies::destroy(condolences->spectral_gen_);
      delete[] condolences->input_spectrum_.data();
      delete[] condolences->input_analyze_.data();
      delete[] condolences->input_buffer_.data();
      delete[] condolences->input_window_.data();
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
      case 3: return smear();
      case 4: return melt();
      case 5: return sensitivity();
      case 6: return decay();
      default: return Parameter::none();
    }
  }

private:
  VESSL_INLINE float get_damping(float in_seconds, const float exp_lin_lerp)
  {
    static constexpr float overlap_size = static_cast<float>(Sympathies::overlap_size);

    // having a shorter decay than the overlap size doesn't make sense
    // and we also want to avoid divide-by-zero.
    in_seconds = vessl::math::max(decay_min_, in_seconds);

    // amplitude needs to decrease by 1 / (decaySeconds * sampleRate()) every sample.
    // eg decaySeconds == 1 -> 1 / sampleRate()
    //    decaySeconds == 0.5 -> 1 / (0.5 * sampleRate), which is twice as fast, equivalent to 2 / sampleRate()
    // since we generate a new buffer every overlapSize samples, we multiply that rate by overlapSize, giving:
    const float damp_lin = overlap_size / (in_seconds * sample_rate_);

    // exponential decay
    const float block_rate = sample_rate_ / overlap_size;
    const float length_in_blocks = in_seconds * block_rate;
    const float damp_exp = 1.0 + vessl::math::log(0.0001f) / (length_in_blocks + 20);
    return vessl::math::lerp(damp_exp, damp_lin, exp_lin_lerp);
  }

  VESSL_INLINE size_t f2t(size_t fbi)
  {
    static constexpr float spectral_band_count = SpectrumSize/2;
    const float tbi_log = spectral_gen_->get_band_index(fbi*band_spacing_) / spectral_band_count;
    const float tbi_lin = vessl::math::interp<vessl::math::easing::expo::out>(0.f, 1.f, tbi_log);
    return static_cast<size_t>(vessl::math::lerp(tbi_log, tbi_lin, spacing_.value) * spectral_band_count);
  } 

  struct
  {
    vessl::analog_p density;
    vessl::analog_p spacing;
    vessl::analog_p spread;
    vessl::analog_p smear;
    vessl::analog_p melt;
    vessl::analog_p response;
    vessl::analog_p decay;
  } params_;
  
  Smoother density_;
  Smoother spacing_;
  Smoother spread_;
  Smoother smear_;
  Smoother melt_;
  Smoother damping_;
  Smoother response_;
  Smoother volume_;
  
  analog_t sample_rate_;
  analog_t band_spacing_;
  analog_t band_first_idx_;
  analog_t band_last_idx_;
  analog_t decay_min_;
  
  SampleArray  input_window_;
  SampleArray  input_buffer_;
  ComplexArray input_spectrum_;
  
  FFT input_fft_;

  Sympathies* spectral_gen_;
  uint32_t    input_buffer_write_idx_;
};