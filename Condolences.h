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
  - expose spread_width as a parameter (integer like density)?
  - try having excite move towards the source spectrum value (pass in actual complex number)
  - revisit volume adjustment formula
*/

template<typename T, uint16_t SpectrumSize, uint8_t Overlap>
class Condolences : public vessl::unit_processor<T>, public vessl::plist<6>
{
public:
  using sample_t     = T;
  using size_t       = vessl::size_t;
  using analog_t     = vessl::analog_t;
  using complex_t    = vessl::transform::complex<sample_t>;
  using Smoother     = vessl::math::easing::smoother<vessl::analog_t>;
  using Parameter    = vessl::parameter;
  using Sympathies   = SpectralSympathies<SpectrumSize, Overlap, false>;
  using SampleArray  = vessl::array<sample_t>;
  using ComplexArray = vessl::array<complex_t>;
  using StringArray  = vessl::array<uint16_t>;
  using Window       = vessl::sample::windows::type;
  using FFT          = vessl::transform::fft<sample_t>;
  using Frequency    = vessl::frequency<analog_t>;

  static constexpr size_t AnalysisSize = Sympathies::overlap_size;
  static constexpr size_t StringCountMax = AnalysisSize/2;
  static constexpr size_t GenerateBlockSize = AnalysisSize;
  static constexpr size_t SpreadWidth = 2;
  // for clamping the param
  static constexpr float DensityMin = 4;
  static constexpr float DensityMax = static_cast<float>(AnalysisSize/2)/(SpreadWidth*2);
  
  Condolences(
    const analog_t sample_rate, 
    sample_t* input_buffer_data, 
    sample_t* input_window_data,
    sample_t* input_analyze_data, 
    complex_t* input_spectrum_data,
    uint16_t*    string_data, 
    Sympathies* spectral_generator
  )
  : sample_rate_(sample_rate)
  , band_spacing_(sample_rate/AnalysisSize)
  , band_first_idx_(1.f + SpreadWidth)
  , band_last_idx_(static_cast<float>(AnalysisSize/2) - SpreadWidth - 1)
  , decay_min_(static_cast<float>(SpectrumSize/2) / sample_rate)
  , input_buffer_(input_buffer_data, AnalysisSize)
  , input_window_(input_window_data, AnalysisSize)
  , input_analyze_(input_analyze_data, AnalysisSize)
  , input_spectrum_(input_spectrum_data, AnalysisSize/2)
  , string_indices_(string_data, StringCountMax)
  , input_fft_(AnalysisSize)
  , spectral_gen_(spectral_generator)
  {
    vessl::sample::windows::render(Window::hann, input_window_);
    input_buffer_.fill(vessl::cast<T>(0.f));
  }

  [[nodiscard]] Parameter density() const { return params_.density("density", 'd'); }
  [[nodiscard]] Parameter spacing() const { return params_.spacing("spacing", 's'); }
  [[nodiscard]] Parameter spread() const { return params_.spread("spread", 'r'); }
  [[nodiscard]] Parameter smear() const { return params_.smear("smear", 'e'); }
  [[nodiscard]] Parameter melt() const { return params_.melt("melt", 'm'); }
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
    //const size_t block_size = in.size();
    
    smear_   = params_.smear.value;
    spread_  = vessl::math::interp<vessl::math::easing::quad::out>(0.f, 1.f, params_.spread.value);
    decay_   = vessl::math::max(decay_min_, params_.decay.value);
    melt_    = params_.melt.value;

    density_ = vessl::math::constrain(params_.density.value, DensityMin, DensityMax);
    spacing_ = params_.spacing.value;
    
    // reduce volume based on combination of decay, spread, and brightness parameters
    // volume_ = vessl::math::interp<vessl::math::easing::expo::out>(1.0f, 0.5f, 
    //     0.2f*params_.decay.value
    //   + 0.2f*params_.spread.value);
    volume_ = 0.25f;
    
    spectral_gen_->spread() = smear_.value;
    spectral_gen_->decay() = vessl::duration_t::from_seconds(decay_.value, sample_rate_);
    spectral_gen_->melt() = melt_.value;
    spectral_gen_->volume() = volume_.value;

    in.copy_to(input_analyze_);


    const size_t string_count = vessl::math::max(static_cast<size_t>(density_.value), 1ull);
    {
      // window the input and output to an analysis buffer
      // because running the fft messes up the input samples.
      input_window_.multiply(input_analyze_, input_analyze_);
      input_fft_.forward(input_analyze_, input_spectrum_);

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
        float damping = 0.9f;
        // main string
        {
          const size_t tbi = spectral_gen_->get_band_index(fbi*band_spacing_); 
          spectral_gen_->excite(tbi, input_spectrum_[fbi], damping);
        }

        // spread strings
        {
          damping *= spread_.value;
          const size_t tbi0 = f2t(fbi-1);
          const size_t tbi1 = f2t(fbi+1);
          spectral_gen_->excite(tbi0, input_spectrum_[fbi-1], damping);
          spectral_gen_->excite(tbi1, input_spectrum_[fbi+1], damping);
        }
        {
          damping *= spread_.value;
          const size_t tbi0 = f2t(fbi-2);
          const size_t tbi1 = f2t(fbi+2);
          spectral_gen_->excite(tbi0, input_spectrum_[fbi-2], damping);
          spectral_gen_->excite(tbi1, input_spectrum_[fbi+2], damping);
        }

        ++si;
        pbi = fbi;
      }
    }

    spectral_gen_->generate(out);
  }

  VESSL_INLINE analog_t get_decay_min() const { return decay_min_; };
  
  VESSL_INLINE typename Sympathies::band_t get_band(analog_t freq_in_hz) const
  {
    return spectral_gen_->get_band(freq_in_hz);
  }
  
  static Condolences* create(vessl::analog_t sample_rate)
  {
    sample_t* input_buffer_data = new sample_t[AnalysisSize];
    sample_t* input_window_data = new sample_t[AnalysisSize];
    sample_t* input_analyze_data = new sample_t[AnalysisSize];
    complex_t* input_spectrum_data = new complex_t[AnalysisSize/2];
    uint16_t*  string_data = new uint16_t[StringCountMax];
    Sympathies* spectral_generator = Sympathies::create(sample_rate);
    return new Condolences(sample_rate,
      input_buffer_data, 
      input_window_data, 
      input_analyze_data,
      input_spectrum_data,
      string_data,
      spectral_generator
      );
  }
  
  static void destroy(const Condolences* condolences)
  {
    if (condolences)
    {
      Sympathies::destroy(condolences->spectral_gen_);
      delete[] condolences->string_indices_.data();
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
      case 3: return smear();
      case 4: return melt();
      case 5: return decay();
      default: return Parameter::none();
    }
  }

private:
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
    vessl::analog_p decay;
    vessl::analog_p feedback;
  } params_;
  
  Smoother density_;
  Smoother spacing_;
  Smoother spread_;
  Smoother smear_;
  Smoother melt_;
  Smoother decay_;
  Smoother volume_;
  
  analog_t sample_rate_;
  analog_t band_spacing_;
  analog_t band_first_idx_;
  analog_t band_last_idx_;
  analog_t decay_min_;
  
  size_t input_buffer_write_;
  SampleArray  input_buffer_;
  SampleArray  input_window_;
  SampleArray  input_analyze_;
  ComplexArray input_spectrum_;
  StringArray  string_indices_;
  
  FFT input_fft_;

  Sympathies* spectral_gen_;
};