#pragma once

#include "vessl/vessl.h"

// lightweight extension that allocates and deallocates the correct amount of memory.
template<typename T, vessl::size_t SpectrumSize, size_t Overlap = 1>
class SpectralGenerator : public vessl::generators::spectral<T, SpectrumSize, Overlap>
{
  using Data = typename vessl::generators::spectral<T, SpectrumSize, Overlap>::data;
  using FrequencyBandType = typename vessl::generators::spectral<T, SpectrumSize, Overlap>::frequency_band;
  using SampleType = typename vessl::generators::spectral<T, SpectrumSize, Overlap>::sample_t;
  using ComplexType = typename vessl::generators::spectral<T, SpectrumSize, Overlap>::complex_t;
  
public:  
  SpectralGenerator(Data& data, vessl::analog_t sample_rate)
    : vessl::generators::spectral<T, SpectrumSize, Overlap>(data, sample_rate)
  {
  }

  static SpectralGenerator* create(vessl::analog_t sample_rate, vessl::sample::windows::type window_type)
  {
    // allocate sample_data first to increase the chances there is a block this big available
    // when SpectrumSize and Overlap are relatively large.
    SampleType* sample_data = new SampleType[SpectrumSize*Overlap*2];
    constexpr vessl::size_t bands_size = SpectrumSize/2;
    FrequencyBandType* bands_data = new FrequencyBandType[bands_size];
    ComplexType* spectrum_data = new ComplexType[bands_size];
    SampleType* window_data = new SampleType[SpectrumSize];
    Data data = {
      vessl::array<FrequencyBandType>(bands_data, bands_size), // bands
      vessl::array<ComplexType>(spectrum_data, bands_size), // spectrum
      vessl::array<SampleType>(sample_data, SpectrumSize*Overlap*2), // signal
      vessl::array<SampleType>(window_data, SpectrumSize), // window
    };
    
    vessl::sample::windows::render(window_type, data.window);
    
    return new SpectralGenerator(data, sample_rate);
  }
  
  static void destroy(SpectralGenerator* generator)
  {
    if (generator)
    {
      delete[] generator->window_.data();
      delete[] generator->signal_[0].data();
      // shouldn't need to delete b because it was allocated along with a in create.
      delete[] generator->spectrum_.data();
      delete[] generator->bands_.data();
    }
    delete generator;
  }
};