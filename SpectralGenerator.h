#pragma once

#include "vessl/vessl.h"

// lightweight extension that allocates and deallocates the correct amount of memory.
template<typename T, vessl::size_t SpectrumSize>
class SpectralGenerator : public vessl::generators::spectral<T, SpectrumSize>
{
  using Data = typename vessl::generators::spectral<T, SpectrumSize>::data;
  using FrequencyBandType = typename vessl::generators::spectral<T, SpectrumSize>::frequency_band;
  using SampleType = typename vessl::generators::spectral<T, SpectrumSize>::sample_t;
  using ComplexType = typename vessl::generators::spectral<T, SpectrumSize>::complex_t;
  
  SpectralGenerator(Data& data, vessl::analog_t sample_rate)
    : vessl::generators::spectral<T, SpectrumSize>(data, sample_rate)
  {
  }
  
public:
  static SpectralGenerator* create(vessl::analog_t sample_rate, vessl::sample::windows::type window_type)
  {
    constexpr vessl::size_t bands_size = SpectrumSize/2;
    FrequencyBandType* frequency_data = new FrequencyBandType[bands_size];
    ComplexType* spectrum_data = new ComplexType[bands_size];
    SampleType* sample_data = new SampleType[SpectrumSize*2];
    SampleType* window_data = new SampleType[SpectrumSize];
    Data data = {
      vessl::array<FrequencyBandType>(frequency_data, bands_size), // frequencies
      vessl::array<ComplexType>(spectrum_data, bands_size), // spectrum
      vessl::array<SampleType>(sample_data, SpectrumSize*2), // signal
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
      delete[] generator->signal_a_.data();
      // shouldn't need to delete b because it was allocated along with a in create.
      delete[] generator->spectrum_.data();
      delete[] generator->frequencies_.data();
    }
    delete generator;
  }
};