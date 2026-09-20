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
class Condolences : public vessl::unit_processor<T>, public vessl::plist<9>
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

  static constexpr size_t AnalysisSize = SpectrumSize/Overlap;

  // for clamping the param
  static constexpr size_t DensityMin = 4;
  static constexpr size_t DensityMax = AnalysisSize/2;

  // frequency in Hz of each string.
  // Density describes how many strings to use.
  // String frequency is then used to exite a specific analysis band.
  // So we pre-compute all string frequencies.
  vessl::analog_t strings[DensityMax];
  
  Condolences(
    const analog_t sample_rate,  
    sample_t* input_window_data,
    sample_t* input_buffer_data,
    sample_t* input_analysis_data,
    complex_t* input_spectrum_data,
    sample_t* output_buffer_data,
    Sympathies* spectral_generator
  )
  : sample_rate_(sample_rate)
  , band_spacing_(sample_rate/AnalysisSize)
  , decay_min_(static_cast<float>(Sympathies::overlap_size) / sample_rate)
  , input_window_(input_window_data, AnalysisSize)
  , input_buffer_(input_buffer_data, AnalysisSize)
  , input_analysis_(input_buffer_data, AnalysisSize)
  , input_spectrum_(input_spectrum_data, AnalysisSize/2)
  , input_fft_(AnalysisSize)
  , spectral_gen_(spectral_generator)
  , output_buffer_(output_buffer_data, Sympathies::overlap_size)
  , input_buffer_write_idx_(0)
  , output_buffer_read_idx_(0)
  {
    params_.response.value = 0.1f;
    input_buffer_.fill(0);
    output_buffer_.fill(0);

    for(size_t i = 0; i < DensityMax; ++i)
    {
      analog_t midi_note = 127*(static_cast<analog_t>(i) / DensityMax);
      strings[i] = vessl::midi_note_to_hertz(midi_note);
    }
  }

  [[nodiscard]] Parameter density() const { return params_.density("density", 'd'); }
  [[nodiscard]] Parameter shift() const { return params_.shift("shift", 'h'); }
  [[nodiscard]] Parameter spacing() const { return params_.spacing("spacing", 's'); }
  [[nodiscard]] Parameter spread() const { return params_.spread("spread", 'r'); }
  [[nodiscard]] Parameter smear() const { return params_.smear("smear", 'e'); }
  [[nodiscard]] Parameter melt() const { return params_.melt("melt", 'm'); }
  [[nodiscard]] Parameter motion() const { return params_.motion("motion", 'o'); }
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
    motion_   = params_.motion.value;
    damping_  = get_damping(params_.decay.value, 0.f);

    density_ = vessl::math::constrain(params_.density.value, static_cast<analog_t>(DensityMin), static_cast<analog_t>(DensityMax));
    shift_   = 1.f + (params_.shift.value < -0.01f ? params_.shift.value*0.75f 
                   : (params_.shift.value > 0.01f ? params_.shift.value*4.f : 0.f));
    spacing_ = params_.spacing.value;
    
    // reduce volume based on combination of decay, spread, and brightness parameters
    // volume_ = vessl::math::interp<vessl::math::easing::expo::out>(1.0f, 0.5f, 
    //     0.2f*params_.decay.value
    //   + 0.2f*params_.spread.value);
    volume_ = 1.0f;
    
    spectral_gen_->smear() = smear_.value;
    spectral_gen_->damping() = damping_.value;
    spectral_gen_->melt() = melt_.value;
    spectral_gen_->motion() = motion_.value;
    spectral_gen_->volume() = volume_.value;

    const size_t string_count = vessl::math::max(static_cast<size_t>(density_.value), 4ull);
    for(size_t i = 0; i < block_size; ++i)
    {
      input_buffer_[input_buffer_write_idx_++] = in[i];

      if (input_buffer_write_idx_ == AnalysisSize)
      {
        input_buffer_.multiply(input_window_, input_analysis_);
        input_fft_.forward(input_analysis_, input_spectrum_);

        // transfer spectrum data from input analysis to spectral_gen
        // by sampling only those frequencies represented by our strings.
        // i.e. comb filter it.
        for(size_t si = 0; si < string_count; ++si)
        {
          const float st = static_cast<float>(si)/(string_count-1);
          const size_t sidx = static_cast<size_t>(st*(DensityMax-1));
          const float shz = strings[sidx];
          const size_t fbi = spectral_gen_->get_band_index(shz);
          const size_t tbi = f2t(fbi, shift_.value);
          float response = response_.value;

          // main string
          excite(tbi, fbi, response);

          // spread strings
          static constexpr size_t ts = SpectrumSize / AnalysisSize;
          {
            response *= spread_.value;
            const size_t fbi0 = fbi-1;
            const size_t fbi1 = fbi+1;
            const size_t tbi0 = tbi - ts; // f2t(fbi-1, shift_.value);
            const size_t tbi1 = tbi + ts; // f2t(fbi+1, shift_.value);
            excite(tbi0, fbi0, response);
            excite(tbi1, fbi1, response);
          }

          {
            response *= spread_.value;
            const size_t fbi0 = fbi-2;
            const size_t fbi1 = fbi+2;
            const size_t tbi0 = tbi - 2*ts; // f2t(fbi-2, shift_.value);
            const size_t tbi1 = tbi + 2*ts; // f2t(fbi+2, shift_.value);
            excite(tbi0, fbi0, response);
            excite(tbi1, fbi1, response);
          }
        }

        /** @todo figure out why this breaks the audio thread */
        // if constexpr(Sympathies::overlap_size < AnalysisSize)
        // {
        //   const size_t count = AnalysisSize - Sympathies::overlap_size;
        //   for(size_t s = 0; s < count; ++s)
        //   {
        //     input_buffer_[s] = input_buffer_[s+Sympathies::overlap_size];
        //   }
        //   input_buffer_write_idx_ = count;
        // }
        // else
        {
          input_buffer_write_idx_ = 0;
        }
      }

      out[i] = output_buffer_[output_buffer_read_idx_++];
      if (output_buffer_read_idx_ == Sympathies::overlap_size)
      {
        spectral_gen_->generate(output_buffer_);
        output_buffer_read_idx_ = 0;
      }
    }
  }

  VESSL_INLINE void excite(const size_t tidx, const size_t fidx, const float response)
  {
    if (!(tidx < 1 || tidx >= SpectrumSize/2 || fidx < 1 || fidx >= AnalysisSize/2))
    {
      spectral_gen_->excite(tidx, input_spectrum_[fidx], response);
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
    sample_t*   output_buffer_data  = new sample_t[Sympathies::overlap_size];
    
    vessl::sample::windows::render(input_window_type, input_window_data, AnalysisSize);
    
    return new Condolences(sample_rate,
      input_window_data,
      input_buffer_data, 
      input_spectrum_data,
      output_buffer_data,
      spectral_generator
      );
  }
  
  static void destroy(const Condolences* condolences)
  {
    if (condolences)
    {
      Sympathies::destroy(condolences->spectral_gen_);
      delete[] condolences->input_spectrum_.data();
      delete[] condolences->output_buffer_.data();
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
      case 1: return shift();
      case 2: return spacing();
      case 3: return spread();
      case 4: return smear();
      case 5: return melt();
      case 6: return sensitivity();
      case 7: return decay();
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

  VESSL_INLINE size_t f2t(size_t fbi, const float thz_scale)
  {
    static constexpr float spectral_band_count = SpectrumSize/2;
    const float thz = (fbi*band_spacing_)*thz_scale;
    const float tbi_log = spectral_gen_->get_band_index(thz) / spectral_band_count;
    const float tbi_lin = vessl::math::interp<vessl::math::easing::expo::out>(0.f, 1.f, tbi_log);
    return static_cast<size_t>(vessl::math::lerp(tbi_log, tbi_lin, spacing_.value) * spectral_band_count);
  } 

  struct
  {
    vessl::analog_p density;
    vessl::analog_p shift;
    vessl::analog_p spacing;
    vessl::analog_p spread;
    vessl::analog_p smear;
    vessl::analog_p melt;
    vessl::analog_p motion;
    vessl::analog_p response;
    vessl::analog_p decay;
  } params_;
  
  Smoother density_;
  Smoother shift_;
  Smoother spacing_;
  Smoother spread_;
  Smoother smear_;
  Smoother melt_;
  Smoother motion_;
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
  SampleArray  input_analysis_;
  ComplexArray input_spectrum_;
  
  FFT input_fft_;

  Sympathies* spectral_gen_;
  SampleArray output_buffer_;
  uint32_t    input_buffer_write_idx_;
  uint32_t    output_buffer_read_idx_;
};