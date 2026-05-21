#pragma once

#include "vessl/vessl.h"

template<typename T = vessl::analog_t>
class Rotator3D : public vessl::unit_processor<vessl::frame::channels<T,3>>
                , protected vessl::plist<10>
{
    using sample_t = T;

public:
    using SampleType = vessl::frame::channels<T,3>;
    using Transform = vessl::transform33<sample_t>;

private:
    using phase_t = vessl::phase_t;
    using phase_p = vessl::phase_p;
    using analog_t = vessl::analog_t;
    using analog_p = vessl::analog_p;
    using q31_t = vessl::q31;
    using q31_p = vessl::param<vessl::q31>;

    Transform rotator;
    phase_t dt; // 1 / sampleRate
    phase_t phaseX, phaseY, phaseZ;

    struct
    {
      // inputs
      analog_p freqInHz;
      analog_p ratioX;
      analog_p ratioY;
      analog_p ratioZ;
      phase_p  modX;
      phase_p  modY;
      phase_p  modZ;

      // outputs
      q31_p rotationX;
      q31_p rotationY;
      q31_p rotationZ;
    } params;

  public:
    explicit Rotator3D(analog_t sampleRate)
      : dt(vessl::cast<phase_t>(1.0f/sampleRate))
      , phaseX(vessl::phase_zero), phaseY(vessl::phase_zero), phaseZ(vessl::phase_zero)
    {

    }

    using param = vessl::parameter;

    [[nodiscard]] VESSL_INLINE param frequency() const { return params.freqInHz({ "frequency", 'f', analog_p::type }); }
    
    // ratios are rotation rates around a particular axis, relative to frequency: [-1,1]
    [[nodiscard]] VESSL_INLINE param ratioX() const { return params.ratioX({"ratio X", 'X', analog_p::type}); }
    [[nodiscard]] VESSL_INLINE param ratioY() const { return params.ratioY({"ratio Y", 'Y', analog_p::type}); }
    [[nodiscard]] VESSL_INLINE param ratioZ() const { return params.ratioZ({"ratio Z", 'Z', analog_p::type}); }
    
    [[nodiscard]] VESSL_INLINE param modX() const   { return params.modX({"mod X", 'x', phase_p::type}); }
    [[nodiscard]] VESSL_INLINE param modY() const   { return params.modY({"mod Y", 'y', phase_p::type}); }
    [[nodiscard]] VESSL_INLINE param modZ() const   { return params.modZ({"mod Z", 'z', phase_p::type}); }

    [[nodiscard]] VESSL_INLINE param rotationX() const { return params.rotationX({"rotation X", 'i', q31_p::type}); }
    [[nodiscard]] VESSL_INLINE param rotationY() const { return params.rotationY({"rotation Y", 'j', q31_p::type}); }
    [[nodiscard]] VESSL_INLINE param rotationZ() const { return params.rotationZ({"rotation Z", 'k', q31_p::type}); }

    [[nodiscard]] const parameter_list& parameters() const override { return *this; }

    VESSL_INLINE void reset() 
    { 
      rotator.set_identity();
      phaseX = 0;
      phaseY = 0;
      phaseZ = 0; 
    }

    VESSL_INLINE SampleType process(const SampleType& in) override
    {
      phase_t rx = phaseX + params.modX.value;
      phase_t ry = phaseY + params.modY.value;
      phase_t rz = phaseZ + params.modZ.value;

      rotator.set_euler(rx, ry, rz);
      
      analog_t freq = params.freqInHz.value;
      analog_t fInc = freq * dt;

      analog_t rxf = params.ratioX.value;
      analog_t ryf = params.ratioY.value;
      analog_t rzf = params.ratioZ.value;

      phaseX += fInc*rxf;
      phaseY += fInc*ryf;
      phaseZ += fInc*rzf;

      params.rotationX.value = vessl::math::sin<q31_t>(rx);
      params.rotationY.value = vessl::math::cos<q31_t>(ry);
      params.rotationZ.value = vessl::math::sin<q31_t>(rz);

      return rotator.process(in);
    }

    VESSL_INLINE void process(vessl::array<SampleType> input, vessl::array<SampleType> output) override
    {
      Transform fromRotator(rotator);

      analog_t freq = params.freqInHz.value;
      analog_t fInc = dt * freq * input.size();

      analog_t rxf = params.ratioX.value;
      analog_t ryf = params.ratioY.value;
      analog_t rzf = params.ratioZ.value;

      phase_t xInc = rxf > 0 ? fInc*rxf : (fInc*rxf) + vessl::phase_360;
      phase_t yInc = ryf > 0 ? fInc*ryf : (fInc*ryf) + vessl::phase_360;
      phase_t zInc = rzf > 0 ? fInc*rzf : (fInc*rzf) + vessl::phase_360;

      phaseX += xInc;
      phaseY += yInc;
      phaseZ += zInc;

      phase_t rx = phaseX + params.modX.value;
      phase_t ry = phaseY + params.modY.value;
      phase_t rz = phaseZ + params.modZ.value;

      rotator.set_euler(rx, ry, rz);

      sample_t* fromData = fromRotator.matrix().data();
      sample_t* toData = rotator.matrix().data();
      sample_t  pct = vessl::cast<sample_t>(1.0f / input.size());

      sample_t rotDeltas[3*3];
      for(int i = 0; i < 9; ++i)
      {
        rotDeltas[i] = (toData[i] - fromData[i]) * pct;
      }

      auto reader = input.make_reader();
      auto writer = output.make_writer();
      while(reader.available())
      {
        SampleType in = reader.read();
        writer << fromRotator.process(in);

        // lerp matrix values towards our target rotation
        for(int i = 0; i < 9; ++i)
        {
          fromData[i] += rotDeltas[i];
        }
      }

      params.rotationX.value = vessl::math::sin<q31_t>(rx);
      params.rotationY.value = vessl::math::cos<q31_t>(ry);
      params.rotationZ.value = vessl::math::sin<q31_t>(rz);

    }

  protected:
    [[nodiscard]] param element_at(vessl::size_t index) const override
    {
      param p[num] = {
        frequency(), 
        ratioX(), ratioY(), ratioZ(),
        modX(), modY(), modZ(),
        rotationX(), rotationY(), rotationZ()
      };
      return p[index];
    }
};