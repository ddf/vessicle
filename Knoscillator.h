#pragma once

#include "vessl/vessl.h"
#include "KnotOscillator.h"
#include "Rotator3D.h"
#include "Noise.hpp"

template<typename T, bool SmoothPQ>
class Knoscillator : public vessl::unit_generator<vessl::sample::frame<T,3>>
  , protected vessl::plist<19>
{
  using sample_t = T;
public:
  using SampleType = vessl::sample::frame<T,3>;
  using KnotOscil = KnotOscillator<sample_t>;
  using KnotType = typename KnotOscil::KnotType;
  
private:
  using SineWave = vessl::sample::waves::bipolar::sine<sample_t>;
  using Rotator = Rotator3D<sample_t>;
  
  using size_t = vessl::size_t;
  using param = vessl::parameter;
  using analog_t = vessl::analog_t;
  using analog_p = vessl::analog_p;
  using phase_t = vessl::phase_t;
  using phase_p = vessl::phase_p;
  using q31_t = vessl::q31;
  using q31_p = vessl::param<q31_t>;
  

public:
  KnotOscil      knoscil;
  Rotator        rotator;
  SineWave       modWave;
  
  phase_t dt;
  phase_t phaseMod;
  
  struct
  {
    // inputs
    analog_p freqInHz;
    analog_p fmRatio;
    analog_p fmIndex;
  } params;

  explicit Knoscillator(float sr)
    : knoscil(sr), rotator(sr)
    , dt(vessl::cast<phase_t>(1.0f/sr))
    , phaseMod(vessl::phase_zero)
  {
    knoscil.knotP() = 2.f;
    knoscil.knotQ() = 1.f;

    rotator.frequency() = 16.f;
    //rotator.ratioY() = (1.f/8.f);
    
    params.fmRatio.value = 2;
    params.fmIndex.value = 0.f;
  }
  
public:
  [[nodiscard]] VESSL_INLINE const KnotOscil& knot() const { return knoscil; }
  
  [[nodiscard]] VESSL_INLINE param knotTypeA() const { return knoscil.knotTypeA(); }
  [[nodiscard]] VESSL_INLINE param knotTypeB() const { return knoscil.knotTypeB(); }
  [[nodiscard]] VESSL_INLINE param knotMorph() const { return knoscil.knotMorph(); }
  [[nodiscard]] VESSL_INLINE param knotP() const { return knoscil.knotP(); }
  [[nodiscard]] VESSL_INLINE param knotQ() const { return knoscil.knotQ(); }
  [[nodiscard]] VESSL_INLINE param knotModP() const { return knoscil.knotModP(); }
  [[nodiscard]] VESSL_INLINE param knotModQ() const { return knoscil.knotModQ(); }
  
  // in Hz
  [[nodiscard]] VESSL_INLINE param frequency() const { return params.freqInHz("frequency", 'f'); }
  [[nodiscard]] VESSL_INLINE param fmRatio() const   { return params.fmRatio("fm ratio", 'R'); }
  [[nodiscard]] VESSL_INLINE param fmIndex() const   { return params.fmIndex("fm index", 'r'); }
  [[nodiscard]] VESSL_INLINE param rotRatioX() const { return rotator.ratioX(); }
  [[nodiscard]] VESSL_INLINE param rotRatioY() const { return rotator.ratioY(); }
  [[nodiscard]] VESSL_INLINE param rotRatioZ() const { return rotator.ratioZ(); }
  [[nodiscard]] VESSL_INLINE param rotModX() const   { return rotator.modX(); }
  [[nodiscard]] VESSL_INLINE param rotModY() const   { return rotator.modY(); }
  [[nodiscard]] VESSL_INLINE param rotModZ() const   { return rotator.modZ(); }
  
  [[nodiscard]] VESSL_INLINE param rotationX() const { return rotator.rotationX(); }
  [[nodiscard]] VESSL_INLINE param rotationY() const { return rotator.rotationY(); }
  [[nodiscard]] VESSL_INLINE param rotationZ() const { return rotator.rotationZ(); }

  [[nodiscard]] const parameter_list& parameters() const override { return *this; }

  VESSL_INLINE void resetRotation()
  { 
    rotator.reset();
  }

  VESSL_INLINE SampleType generate() override
  { 
    analog_t freq = params.freqInHz.value;
    // phase modulate in sync with the current frequency
    analog_t fmRatio = params.fmRatio.value;
    sample_t fmIndex = vessl::cast<sample_t>(params.fmIndex.value);
    phase_t  mInc = (freq*fmRatio) * dt;
    
    phase_t fm = vessl::cast<phase_t>(modWave.evaluate(phaseMod)*fmIndex);
    phaseMod += mInc;
    
    knoscil.frequency() = freq;
    knoscil.phaseMod()  = fm;

    SampleType coord = knoscil.template generate<SmoothPQ>();
    coord = rotator.process(coord);
  
    return coord;
  }

  VESSL_INLINE void generate(vessl::array<SampleType> dest)
  { 
    analog_t freq = params.freqInHz.value;
    // phase modulate in sync with the current frequency
    analog_t fmRatio = params.fmRatio.value;
    sample_t fmIndex = vessl::cast<sample_t>(params.fmIndex.value);
    phase_t  mInc = (freq*fmRatio) * dt;
    
    knoscil.frequency() = freq;

    // generate knot, then rotate it.
    auto writer = dest.make_writer();
    while(writer.available())
    {
      phase_t fm = vessl::cast<phase_t>(modWave.evaluate(phaseMod)*fmIndex);
      phaseMod += mInc;

      knoscil.phaseMod() = fm;
      writer << knoscil.template generate<SmoothPQ>();
    }
    
    rotator.process(dest, dest);
  }
  
  static Knoscillator* create(float sampleRate)
  {
    return new Knoscillator(sampleRate);
  }
  
  static void destroy(const Knoscillator* knoscillator)
  {
    delete knoscillator;
  }
  
protected:
  VESSL_INLINE param element_at(vessl::size_t index) const override
  {
    param p[plist::num] = {
      knotTypeA(), knotTypeB(), knotMorph(),
      knotP(), knotQ(), knotModP(), knotModQ(),
      frequency(), fmRatio(), fmIndex(), 
      rotRatioX(), rotRatioY(), rotRatioZ(),
      rotModX(), rotModY(), rotModZ(),
      rotationX(), rotationY(), rotationZ()
    };
    return p[index];
  }
};
