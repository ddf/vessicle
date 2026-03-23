#pragma once

#include "KnotOscillator.h"
#include "Noise.hpp"
#include "vessl/vessl.h"

template<typename T = vessl::analog_t>
class Knoscillator : public vessl::unitGenerator<vessl::frame::channels<T, 2>>
  , protected vessl::plist<22>
{
  using sample_t = T;
public:
  using SampleType = vessl::frame::channels<sample_t, 2>;
  using KnotType = typename KnotOscillator<sample_t>::KnotType;
  
private:
  using SineWave = vessl::waves::sine<sample_t>;
  using KnotOscil = KnotOscillator<sample_t>;
  using Smoother = vessl::smoother<>;
  using Transform = vessl::transform33<sample_t>;
  
  using size_t = vessl::size_t;
  using param = vessl::parameter;
  using analog_t = vessl::analog_t;
  using analog_p = vessl::analog_p;
  using coord_t = typename KnotOscil::coord_t;
  using phase_t = vessl::phase_t;
  using phase_p = vessl::phase_p;
  using q31_t = vessl::q31;
  
  static constexpr size_t noiseDim = 128;
  static constexpr float  noiseStep = 4.0f / noiseDim;
  static constexpr analog_t zoomFar = 60.0f * KnotOscil::KNOT_SCALE;
  static constexpr analog_t zoomNear = 6.0f * KnotOscil::KNOT_SCALE;
  
  using NoiseTable = vessl::wavetable<float, noiseDim*noiseDim>;

  KnotOscil      knoscil;
  Transform      rotator;
  Smoother       zoom;
  SineWave       modWave;
  
  phase_t dt;
  phase_t phaseMod;
  phase_t phaseS;
  phase_t rotateX;
  phase_t rotateY;
  phase_t rotateZ;
  T       projection;
  
  struct
  {
    // inputs
    analog_p freqInHz;
    analog_p fmRatio;
    analog_p fmIndex;
    phase_p rotRatioX;
    phase_p rotRatioY;
    phase_p rotRatioZ;
    phase_p rotModX;
    phase_p rotModY;
    phase_p rotModZ;
    phase_p zoom;
    analog_p squiggleAmt;
    analog_p noiseAmt;
    
    // outputs
    analog_p rotationX;
    analog_p rotationY;
    analog_p rotationZ;
  } params;
  
  NoiseTable noiseTable;

  explicit Knoscillator(float sr)
    : knoscil(sr)
    , zoom(0.9f, zoomNear)
    , dt(vessl::cast<phase_t>(1.0f/sr)), phaseMod(vessl::PHASE_ZERO), phaseS(vessl::PHASE_ZERO)
    , rotateX(vessl::PHASE_ZERO), rotateY(vessl::PHASE_ZERO), rotateZ(vessl::PHASE_ZERO)
  {
    knoscil.knotP() = 2.f;
    knoscil.knotQ() = 1.f;
    
    params.fmRatio.value = 2;
    params.fmIndex.value = 0.f;
    params.zoom.value = vessl::PHASE_360;
    params.rotRatioY.value = vessl::cast<phase_t>(1.0f/8.0f);
    
    // for (size_t x = 0; x < noiseDim; ++x)
    // {
    //   for (size_t y = 0; y < noiseDim; ++y)
    //   {
    //     size_t i = x * noiseDim + y;
    //     noiseTable.set(i, vessicle::perlin2d(x*noiseStep, y*noiseStep, 1, 4) * 2 - 1);
    //   }
    // }
  }
  
public:
  const KnotOscil& knot() const { return knoscil; }
  
  param knotTypeA() const { return knoscil.knotTypeA(); }
  param knotTypeB() const { return knoscil.knotTypeB(); }
  param knotMorph() const { return knoscil.knotMorph(); }
  param knotP() const { return knoscil.knotP(); }
  param knotQ() const { return knoscil.knotQ(); }
  param knotModP() const { return knoscil.knotModP(); }
  param knotModQ() const { return knoscil.knotModQ(); }
  
  // in Hz
  param frequency() const { return params.freqInHz({ "frequency", 'f', analog_p::type }); }
  param fmRatio() const   { return params.fmRatio({"fm ratio", 'R', analog_p::type }); }
  param fmIndex() const   { return params.fmIndex({"fm index", 'r', analog_p::type}); }
  param rotRatioX() const { return params.rotRatioX({"rotation ratio X", 'X', analog_p::type}); }
  param rotRatioY() const { return params.rotRatioY({"rotation ratio Y", 'Y', analog_p::type}); }
  param rotRatioZ() const { return params.rotRatioZ({"rotation ratio Z", 'Z', analog_p::type}); }
  param rotModX() const   { return params.rotModX({"rotation mod X", 'x', phase_p::type }); }
  param rotModY() const   { return params.rotModY({"rotation mod Y", 'y', phase_p::type }); }
  param rotModZ() const   { return params.rotModZ({"rotation mod Z", 'z', phase_p::type}); }
  param cameraZoom() const{ return params.zoom({"camera zoom", 'C', analog_p::type}); }
  param squiggle() const  { return params.squiggleAmt({"squiggle amount", 'S', analog_p::type}); }
  param noise() const     { return params.noiseAmt({"noise amount", 'N', analog_p::type}); }
  
  param rotationX() const { return params.rotationX({"rotation X", 'i', analog_p::type}); }
  param rotationY() const { return params.rotationY({"rotation Y", 'j', analog_p::type}); }
  param rotationZ() const { return params.rotationZ({"rotation Z", 'k', analog_p::type}); }

  [[nodiscard]] const parameters& getParameters() const override { return *this; }

  T getProjection() const { return projection; }

  SampleType generate() override
  {
    SampleType out;
    zoom = vessl::easing::lerpp(zoomFar, zoomNear, params.zoom.value);
    //zoom = zoomFar + (zoomNear - zoomFar)*params.zoom.value;

    analog_t sVol = params.squiggleAmt.value * 0.25f;

    phase_t rxm = params.rotModX.value;
    phase_t rym = params.rotModY.value;
    phase_t rzm = params.rotModZ.value;

    analog_t nVol = params.noiseAmt.value * 0.5f;
    
    analog_t freq = params.freqInHz.value;
    // phase modulate in sync with the current frequency
    analog_t fmRatio = params.fmRatio.value;
    sample_t fmIndex = vessl::cast<sample_t>(params.fmIndex.value);
    phase_t mInc = (freq*fmRatio) * dt;
    
    phase_t fm = vessl::cast<phase_t>(modWave.evaluate(phaseMod)*fmIndex);
    phaseMod += mInc;
    
    knoscil.frequency() = freq;
    knoscil.phaseMod()  = fm;

    coord_t coord = knoscil.template generate<false>();
    
    rotator.setEuler(rotateX + rxm, rotateY + rym, rotateZ + rzm);
    coord = rotator.process(coord);
    
    // phase_t st = phaseS + fm;
    // analog_t nz = nVol * noise(coord.x, coord.y);
    // coord.x += vessl::math::cos<analog_t>(st)*sVol + coord.x * nz;
    // coord.y += vessl::math::sin<analog_t>(st)*sVol + coord.y * nz;
    // coord.z += coord.z * nz;

    analog_t zm = zoom.value;
    analog_t cz = vessl::cast<analog_t>(coord.z);
    projection = vessl::cast<sample_t>(1.0f / (cz+zm));
    out.left()  = (coord.x * projection);
    out.right() = (coord.y * projection);
    
    q31_t fInc = q31_t((int64_t)(freq * dt));
    analog_t knotP = knoscil.knotP().readAnalog();
    analog_t knotQ = knoscil.knotQ().readAnalog();
    phase_t sInc = static_cast<phase_t>(fInc * 4 * (knotP + knotQ));
    phaseS  = phaseS + sInc;
    
    // float knotP = knoscil.knotP().readAnalog();
    // float knotQ = knoscil.knotQ().readAnalog();

    q31_t rxf = vessl::cast<q31_t>(params.rotRatioX.value);
    q31_t ryf = vessl::cast<q31_t>(params.rotRatioY.value);
    q31_t rzf = vessl::cast<q31_t>(params.rotRatioZ.value);
    
    // phase_t sInc  = static_cast<phase_t>(fInc * 4 * (knotP + knotQ));
    // phaseS  = phaseS + static_cast<phase_t>(sInc);
    
    q31_t rInc = fInc;
    rotateX += vessl::cast<phase_t>(rInc*rxf);
    rotateY += vessl::cast<phase_t>(rInc*ryf);
    rotateZ += vessl::cast<phase_t>(rInc*rzf);
    
    params.rotationX.value = vessl::math::sin<analog_t>(rotateX + rxm);
    params.rotationY.value = vessl::math::cos<analog_t>(rotateY + rym);
    params.rotationZ.value = vessl::math::sin<analog_t>(rotateZ + rzm);
  
    return out;
  }

  void generate(vessl::array<SampleType> dest)
  {
    // float sVol = params.squiggleAmt.value * 0.25f;
    // float nVol = params.noiseAmt.value * 0.5f;
    
    analog_t freq = params.freqInHz.value;
    // phase modulate in sync with the current frequency
    analog_t fmRatio = params.fmRatio.value;
    sample_t fmIndex = vessl::cast<sample_t>(params.fmIndex.value);
    phase_t mInc = (freq*fmRatio) * dt;
    
    knoscil.frequency() = freq;

    phase_t rxm = params.rotModX.value;
    phase_t rym = params.rotModY.value;
    phase_t rzm = params.rotModZ.value;
    rotator.setEuler(rotateX + rxm, rotateY + rym, rotateZ + rzm);
    T zm = vessl::cast<T>(zoomNear);

    SampleType out;
    auto writer = dest.getWriter();
    while(writer.available())
    {
      phase_t fm = vessl::cast<phase_t>(modWave.evaluate(phaseMod)*fmIndex);
      phaseMod += mInc;

      knoscil.phaseMod() = fm;
      coord_t coord = knoscil.template generate<false>();
      coord = rotator.process(coord);
      
      // phase_t st = phaseS + fm;
      // float nz = nVol * noise(coord.x, coord.y);
      // coord.x += vessl::math::cosz<float>(st)*sVol + coord.x * nz;
      // coord.y += vessl::math::sinz<float>(st)*sVol + coord.y * nz;
      // coord.z += coord.z * nz;

      //analog_t zm = zoomNear; // zoom.value;
      //analog_t cz = vessl::cast<analog_t>(coord.z);
      projection =  zm + coord.z;
      out.left()  = coord.x / projection;
      out.right() = coord.y / projection;

      writer << out;
    }

    // float knotP = knoscil.knotP().readAnalog();
    // float knotQ = knoscil.knotQ().readAnalog();
    q31_t fInc = q31_t((int64_t)(freq * dt));
    q31_t rxf = vessl::cast<q31_t>(params.rotRatioX.value);
    q31_t ryf = vessl::cast<q31_t>(params.rotRatioY.value);
    q31_t rzf = vessl::cast<q31_t>(params.rotRatioZ.value);
    
    // phase_t sInc  = static_cast<phase_t>(fInc * 4 * (knotP + knotQ));
    // phaseS  = phaseS + static_cast<phase_t>(sInc);
    
    q31_t rInc = fInc;
    rotateX += vessl::cast<phase_t>(rInc*rxf);
    rotateY += vessl::cast<phase_t>(rInc*ryf);
    rotateZ += vessl::cast<phase_t>(rInc*rzf);
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
  param elementAt(vessl::size_t index) const override
  {
    param p[plsz] = {
      knotTypeA(), knotTypeB(), knotMorph(), knotP(), knotQ(), knotModP(), knotModQ(),
      frequency(), fmRatio(), fmIndex(), rotRatioX(), rotRatioY(), rotRatioZ(),
      rotModX(), rotModY(), rotModZ(), cameraZoom(), squiggle(), noise(),
      rotationX(), rotationY(), rotationZ()
    };
    return p[index];
  }
  
private:
  [[nodiscard]] float noise(float x, float y) const
  {
    size_t nx = static_cast<size_t>(vessl::math::abs(x) / noiseStep) % noiseDim;
    size_t ny = static_cast<size_t>(vessl::math::abs(y) / noiseStep) % noiseDim;
    size_t ni = nx * noiseDim + ny;
    return noiseTable.get(ni);
  }
};
