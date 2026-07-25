#pragma once

#include "vessl/vessl.h"
#include "Window.h"

// lightweight extension that allocates and deallocates the correct amount of memory.
template<typename T, vessl::size_t SpectrumSize, vessl::size_t Overlap = 2>
class SpectralGenerator : public vessl::generators::spectral<T, SpectrumSize, Overlap>
{
  using Data = typename vessl::generators::spectral<T, SpectrumSize, Overlap>::data;
  using FrequencyBandType = typename vessl::generators::spectral<T, SpectrumSize, Overlap>::frequency_band;
  using SampleType = typename vessl::generators::spectral<T, SpectrumSize, Overlap>::sample_t;
  using ComplexType = typename vessl::generators::spectral<T, SpectrumSize, Overlap>::complex_t;
  
  SpectralGenerator(Data& data, vessl::analog_t sample_rate)
    : vessl::generators::spectral<T, SpectrumSize, Overlap>(data, sample_rate)
  {
  }
  
public:
  static SpectralGenerator* create(vessl::analog_t sample_rate, vessl::sample::windows::type window_type)
  {
    constexpr vessl::size_t bands = SpectrumSize/2;
    Data data = {
      vessl::array<FrequencyBandType>(new FrequencyBandType[bands], bands),
      vessl::array<ComplexType>(new ComplexType[SpectrumSize], SpectrumSize),
      vessl::array<SampleType>(new SampleType[SpectrumSize], SpectrumSize),
      vessl::array<SampleType>(new SampleType[SpectrumSize], SpectrumSize),
      vessl::sample::ring_buffer<SampleType>(new SampleType[SpectrumSize*8], SpectrumSize*8)
    };
    
    //Window::triangular(data.window.data(), data.window.size());
    vessl::sample::windows::render(window_type, data.window);
    
    return new SpectralGenerator(data, sample_rate);
  }
  
  static void destroy(SpectralGenerator* generator)
  {
    if (generator)
    {
      delete[] generator->buffer_.data();
      delete[] generator->window_.data();
      delete[] generator->signal_.data();
      delete[] generator->spectrum_.data();
      delete[] generator->frequencies_.data();
    }
    delete generator;
  }
};