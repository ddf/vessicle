#pragma once

#include "KnotOscillator.h"
#include "Rotator3D.h"
#include "Noise.hpp"
#include "vessl/vessl.h"

template<typename T = vessl::analog_t, vessl::size_t BLOCK_SIZE = 1, bool smooth_pq = true>
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
  using Rotator = Rotator3D<sample_t>;
  
  using size_t = vessl::size_t;
  using param = vessl::parameter;
  using analog_t = vessl::analog_t;
  using analog_p = vessl::analog_p;
  using coord_t = typename KnotOscil::coord_t;
  using phase_t = vessl::phase_t;
  using phase_p = vessl::phase_p;
  using q31_t = vessl::q31;
  using q31_p = vessl::param<q31_t>;
  
  static constexpr size_t noiseDim = 128;
  static constexpr float  noiseStep = 4.0f / noiseDim;
  static constexpr analog_t zoomFar = 60.0f * KnotOscil::KNOT_SCALE;
  static constexpr analog_t zoomNear = 6.0f * KnotOscil::KNOT_SCALE;
  
  using NoiseTable = vessl::wavetable<float, noiseDim*noiseDim>;

public:
  KnotOscil      knoscil;
  Rotator        rotator;
  Smoother       zoom;
  SineWave       modWave;
  
  phase_t dt;
  phase_t phaseMod;
  phase_t phaseS;
  T       projection;
  
  struct
  {
    // inputs
    analog_p freqInHz;
    analog_p fmRatio;
    analog_p fmIndex;
    phase_p zoom;
    analog_p squiggleAmt;
    analog_p noiseAmt;
  } params;

  coord_t blockBuffer[BLOCK_SIZE];
  
  NoiseTable noiseTable;

  explicit Knoscillator(float sr)
    : knoscil(sr), rotator(sr)
    , zoom(0.9f, zoomNear)
    , dt(vessl::cast<phase_t>(1.0f/sr)), phaseMod(vessl::PHASE_ZERO), phaseS(vessl::PHASE_ZERO)
  {
    knoscil.knotP() = 2.f;
    knoscil.knotQ() = 1.f;

    rotator.frequency() = 16.f;
    //rotator.ratioY() = (1.f/8.f);
    
    params.fmRatio.value = 2;
    params.fmIndex.value = 0.f;
    params.zoom.value = vessl::PHASE_360;
    
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
  coord_t knotCoord;
  coord_t knotCoordRotated;

  [[nodiscard]] VESSL_INLINE const KnotOscil& knot() const { return knoscil; }
  
  [[nodiscard]] VESSL_INLINE param knotTypeA() const { return knoscil.knotTypeA(); }
  [[nodiscard]] VESSL_INLINE param knotTypeB() const { return knoscil.knotTypeB(); }
  [[nodiscard]] VESSL_INLINE param knotMorph() const { return knoscil.knotMorph(); }
  [[nodiscard]] VESSL_INLINE param knotP() const { return knoscil.knotP(); }
  [[nodiscard]] VESSL_INLINE param knotQ() const { return knoscil.knotQ(); }
  [[nodiscard]] VESSL_INLINE param knotModP() const { return knoscil.knotModP(); }
  [[nodiscard]] VESSL_INLINE param knotModQ() const { return knoscil.knotModQ(); }
  
  // in Hz
  [[nodiscard]] VESSL_INLINE param frequency() const { return params.freqInHz({ "frequency", 'f', analog_p::type }); }
  [[nodiscard]] VESSL_INLINE param fmRatio() const   { return params.fmRatio({"fm ratio", 'R', analog_p::type }); }
  [[nodiscard]] VESSL_INLINE param fmIndex() const   { return params.fmIndex({"fm index", 'r', analog_p::type}); }
  [[nodiscard]] VESSL_INLINE param rotRatioX() const { return rotator.ratioX(); }
  [[nodiscard]] VESSL_INLINE param rotRatioY() const { return rotator.ratioY(); }
  [[nodiscard]] VESSL_INLINE param rotRatioZ() const { return rotator.ratioZ(); }
  [[nodiscard]] VESSL_INLINE param rotModX() const   { return rotator.modX(); }
  [[nodiscard]] VESSL_INLINE param rotModY() const   { return rotator.modY(); }
  [[nodiscard]] VESSL_INLINE param rotModZ() const   { return rotator.modZ(); }
  [[nodiscard]] VESSL_INLINE param cameraZoom() const{ return params.zoom({"camera zoom", 'C', analog_p::type}); }
  [[nodiscard]] VESSL_INLINE param squiggle() const  { return params.squiggleAmt({"squiggle amount", 'S', analog_p::type}); }
  [[nodiscard]] VESSL_INLINE param noise() const     { return params.noiseAmt({"noise amount", 'N', analog_p::type}); }
  
  [[nodiscard]] VESSL_INLINE param rotationX() const { return rotator.rotationX(); }
  [[nodiscard]] VESSL_INLINE param rotationY() const { return rotator.rotationY(); }
  [[nodiscard]] VESSL_INLINE param rotationZ() const { return rotator.rotationZ(); }

  [[nodiscard]] const parameters& getParameters() const override { return *this; }

  [[nodiscard]] VESSL_INLINE T getProjection() const { return projection; }

  VESSL_INLINE void resetRotation()
  { 
    rotator.reset();
  }

  VESSL_INLINE SampleType generate() override
  {
    SampleType out;
    zoom = vessl::easing::lerpp(zoomFar, zoomNear, params.zoom.value);
    //zoom = zoomFar + (zoomNear - zoomFar)*params.zoom.value;

    //analog_t sVol = params.squiggleAmt.value * 0.25f;
    //analog_t nVol = params.noiseAmt.value * 0.5f;
    
    analog_t freq = params.freqInHz.value;
    // phase modulate in sync with the current frequency
    analog_t fmRatio = params.fmRatio.value;
    sample_t fmIndex = vessl::cast<sample_t>(params.fmIndex.value);
    phase_t mInc = (freq*fmRatio) * dt;
    
    phase_t fm = vessl::cast<phase_t>(modWave.evaluate(phaseMod)*fmIndex);
    phaseMod += mInc;
    
    knoscil.frequency() = freq;
    knoscil.phaseMod()  = fm;

    coord_t coord = knoscil.template generate<smooth_pq>();
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
    
    // phase_t sInc  = static_cast<phase_t>(fInc * 4 * (knotP + knotQ));
    // phaseS  = phaseS + static_cast<phase_t>(sInc);
  
    return out;
  }

  VESSL_INLINE void generate(vessl::array<SampleType> dest)
  {
    // float sVol = params.squiggleAmt.value * 0.25f;
    // float nVol = params.noiseAmt.value * 0.5f;
    
    analog_t freq = params.freqInHz.value;
    // phase modulate in sync with the current frequency
    analog_t fmRatio = params.fmRatio.value;
    sample_t fmIndex = vessl::cast<sample_t>(params.fmIndex.value);
    phase_t mInc = (freq*fmRatio) * dt;
    
    knoscil.frequency() = freq;

    //rotator.setEuler(0, vessl::PHASE_90, 0);
    T zm = vessl::cast<T>(zoomNear);

    // generate knot, then rotate it.
    vessl::array<coord_t> block(blockBuffer, BLOCK_SIZE);
    auto writer = block.getWriter();
    while(writer.available())
    {
      phase_t fm = vessl::cast<phase_t>(modWave.evaluate(phaseMod)*fmIndex);
      phaseMod += mInc;

      knoscil.phaseMod() = fm;
      writer << knoscil.template generate<smooth_pq>();
    }
    
    rotator.process(block, block);

    // project to 2D
    coord_t coord;
    SampleType out;
    auto input = block.getReader();
    auto output = dest.getWriter();
    while(output.available())
    {
      coord = input.read();

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

      output << out;
    }

    // float knotP = knoscil.knotP().readAnalog();
    // float knotQ = knoscil.knotQ().readAnalog();
    // q31_t fInc = q31_t((int64_t)(dest.getSize() * dt));
    // phase_t sInc  = static_cast<phase_t>(fInc * 4 * (knotP + knotQ));
    // phaseS  = phaseS + static_cast<phase_t>(sInc);
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
  VESSL_INLINE param elementAt(vessl::size_t index) const override
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
  [[nodiscard]] VESSL_INLINE float noise(float x, float y) const
  {
    size_t nx = static_cast<size_t>(vessl::math::abs(x) / noiseStep) % noiseDim;
    size_t ny = static_cast<size_t>(vessl::math::abs(y) / noiseStep) % noiseDim;
    size_t ni = nx * noiseDim + ny;
    return noiseTable.get(ni);
  }
};
