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

  static constexpr size_t overlap_size = Sympathies::overlap_size;
  static constexpr size_t generate_block_size = overlap_size;
  static constexpr size_t spread_width = 2;
  // for clamping the param
  static constexpr float density_min = 4;
  static constexpr float density_max = SpectrumSize/spread_width;
  
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
  , band_first_idx_(1.f + spread_width)
  , band_last_idx_(static_cast<float>(SpectrumSize/2) - spread_width - 1)
  , decay_min_(static_cast<float>(SpectrumSize) * 0.5f / sample_rate)
  , input_buffer_write_(SpectrumSize - overlap_size)
  , input_buffer_(input_buffer_data, SpectrumSize)
  , input_window_(input_window_data, SpectrumSize)
  , input_analyze_(input_analyze_data, SpectrumSize)
  , input_spectrum_(input_spectrum_data, SpectrumSize/2)
  , string_indices_(string_data, SpectrumSize/2)
  , input_transform_(SpectrumSize)
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
    const size_t block_size = in.size();
    
    smear_   = params_.smear.value;
    spread_  = vessl::math::interp<vessl::math::easing::quad::out>(0.f, 1.f, params_.spread.value);
    decay_   = vessl::math::max(decay_min_, params_.decay.value);
    melt_    = params_.melt.value;

    density_ = vessl::math::constrain(params_.density.value, density_min, density_max);
    spacing_ = params_.spacing.value;
    
    const size_t string_count = vessl::math::max(static_cast<size_t>(density_.value + 0.5f), 1ull);
    //if (string_count != string_count_ || vessl::math::abs(ps - spacing_.value) > 0.01f)
    {
      uint16_t pbi = 0;
      for (size_t si = 0; si < string_count; ++si)
      {
        const float st = static_cast<float>(si)/string_count;
        const float lin_bi = vessl::math::lerp(band_first_idx_, band_last_idx_, st);
        const float exp_bi = vessl::math::interp<vessl::math::easing::expo::in>(band_first_idx_, band_last_idx_, st);
        const float cbi = vessl::math::lerp(exp_bi, lin_bi, spacing_.value);
        pbi = string_indices_[si] = cbi > pbi ? cbi : pbi+1;
      }
    }
    string_count_ = string_count;
    
    // reduce volume based on combination of decay, spread, and brightness parameters
    // volume_ = vessl::math::interp<vessl::math::easing::expo::out>(1.0f, 0.5f, 
    //     0.2f*params_.decay.value
    //   + 0.2f*params_.spread.value);
    volume_ = 0.5f;
    
    spectral_gen_->spread() = smear_.value;
    spectral_gen_->decay() = vessl::duration_t::from_seconds(decay_.value, sample_rate_);
    spectral_gen_->melt() = melt_.value;
    spectral_gen_->volume() = volume_.value;

    SampleArray buffer(input_buffer_.data() + input_buffer_write_, block_size);
    in.copy_to(buffer);
    input_buffer_write_ += block_size;

    // input buffer will be equal to overlap size
    // and we start writing overlap size from the end of the buffer.
    // so every block we can update our spectrum.
    //if (input_buffer_write_ == SpectrumSize)
    {
      // window the input and output to an analysis buffer
      // because running the fft messes up the input samples.
      input_window_.multiply(input_buffer_, input_analyze_);
      input_transform_.forward(input_analyze_, input_spectrum_);
      
      // copy the back section of the array to the front section
      // continue recording input from overlap_size before the end.
      // doing this means we can update the spectral data for sound generation every overlap.
      input_buffer_write_ = SpectrumSize - overlap_size;
      SampleArray input_buffer_back(input_buffer_.data() + overlap_size, input_buffer_write_);
      input_buffer_back.copy_to(input_buffer_);

      // transfer spectrum data from input analysis to spectral_gen
      // by sampling only those frequencies represented by our strings.
      // i.e. comb filter it.
      size_t si = 0;
      while(si < string_count)
      {
        const size_t bi = string_indices_[si];
        //if (bi > 0 && bi < iss)
        {
          float damping = 0.1f;
          spectral_gen_->excite(bi, input_spectrum_[bi], damping);

          //float in_mag = input.magnitude() * mag;
          //vessl::phase_t in_phase = input.phase();
          //spectral_gen_->excite(bi, in_mag, in_phase);
          
          //size_t si = 0;
          //while(si++ < spread_width)
          // hand unroll to hopefully speed this up
          {
            damping *= spread_.value;
            spectral_gen_->excite(bi-1, input_spectrum_[bi-1], damping);
            spectral_gen_->excite(bi+1, input_spectrum_[bi+1], damping);
          }
          {
            damping *= spread_.value;
            spectral_gen_->excite(bi-2, input_spectrum_[bi-2], damping);
            spectral_gen_->excite(bi+2, input_spectrum_[bi+2], damping);
          }
        }
        ++si;
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
    sample_t* input_buffer_data = new sample_t[SpectrumSize];
    sample_t* input_window_data = new sample_t[SpectrumSize];
    sample_t* input_analyze_data = new sample_t[SpectrumSize];
    complex_t* input_spectrum_data = new complex_t[SpectrumSize/2];
    uint16_t*  string_data = new uint16_t[SpectrumSize/2];
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
  analog_t band_first_idx_;
  analog_t band_last_idx_;
  analog_t decay_min_;
  
  size_t input_buffer_write_;
  SampleArray  input_buffer_;
  SampleArray  input_window_;
  SampleArray  input_analyze_;
  ComplexArray input_spectrum_;
  StringArray  string_indices_;
  size_t       string_count_ = 0;
  size_t       string_process_index_ = 0;
  
  FFT input_transform_;

  Sympathies* spectral_gen_;
};