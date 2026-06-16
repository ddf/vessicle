#pragma once

#include "vessl/vessl.h"

template<typename T = vessl::analog_t>
class KnotOscillator : public vessl::unit_generator<vessl::sample::frame<T,3>>
                     , protected vessl::plist<9>
{
public:
  enum class KnotType : uint8_t
  {
    TFOIL = 0,
    LISSA = 1,
    TORUS = 2,

    COUNT = 3 // note: update interp method if more knots are added
  };

  using parameter = vessl::parameter;
  using sample_t  = T;
  using coord_t   = vessl::sample::frame<T,3>;
  using phase_t   = vessl::phase_t;
  using phasew_t  = vessl::size_t;
  using analog_t  = vessl::analog_t;
  using digital_p = vessl::digital_p;
  using analog_p  = vessl::analog_p;
  using knot_p    = vessl::param<KnotType>;
  using phase_p   = vessl::phase_p;
  
  static constexpr int KNOT_TYPE_COUNT = static_cast<int>(KnotType::COUNT);
  static constexpr analog_t KNOT_SCALE = (1.f / 12.f);
  static constexpr sample_t TORUS_SCALE = vessl::cast<sample_t>(0.2f);
  
private:
  struct
  {
    knot_p knotTypeA;
    knot_p knotTypeB;
    phase_p  knotMorph;
    analog_p knotP;
    analog_p knotQ;
    analog_p frequency;
    phase_p  phaseMod;
    analog_p knotModP;
    analog_p knotModQ;
  } params;

  // x3 and y2 are used as phase_t values, but in order for morphing to work,
  // we need to be able to morph from 0 to 3*PI radians, which isn't possible with phase_t.
  // morphing only from 0 to PI is not sufficient to line up the different knots correctly.
  // all other coefficients use type T so that we can support fixed-point types.
  // for this same reason, all T coefficients are in the range [-1,1]
  sample_t x1[KNOT_TYPE_COUNT];
  sample_t x2[KNOT_TYPE_COUNT];
  phasew_t x3[KNOT_TYPE_COUNT];
  sample_t y1[KNOT_TYPE_COUNT];
  phasew_t y2[KNOT_TYPE_COUNT];
  sample_t y3[KNOT_TYPE_COUNT];
  sample_t z1[KNOT_TYPE_COUNT];
  sample_t z2[KNOT_TYPE_COUNT];

  // cached coefficient values that don't depend on phase.
  // updated only when knot types change or morph changes.
  sample_t cx1;
  phase_t  cx3;
  sample_t cy1;
  phase_t  cy2;
  sample_t cz1;
  sample_t cz2;

  KnotType knotA = KnotType::COUNT;
  KnotType knotB = KnotType::COUNT;
  phase_t  knotM = 0;
  
  phase_t phaseP;
  phase_t phaseQ;
  phase_t phaseZ;
  phase_t dt;
  float   sr_;
  coord_t xyz_;

public:
  explicit KnotOscillator(float sampleRate)
    : params()
    , phaseP(vessl::phase_zero), phaseQ(vessl::phase_zero), phaseZ(vessl::phase_zero)
    , dt(vessl::cast<phase_t>(1.0f/sampleRate)), sr_(sampleRate)
  {
    params.knotP.value = 1;
    params.knotQ.value = 1;
    params.knotTypeA.value = KnotType::TFOIL;
    params.knotTypeB.value = KnotType::LISSA;
    params.frequency.value = 1;
    
    static constexpr int TFOIL = static_cast<int>(KnotType::TFOIL);
    x1[TFOIL] = vessl::cast<sample_t>(1.f * KNOT_SCALE);
    x2[TFOIL] = vessl::cast<sample_t>(2.f * KNOT_SCALE);
    x3[TFOIL] = 3 * vessl::cast<phasew_t>(vessl::phase_90); // 3*PI/2;
    y1[TFOIL] = vessl::cast<sample_t>(1.f * KNOT_SCALE);
    y2[TFOIL] = vessl::phase_zero;
    y3[TFOIL] = vessl::cast<sample_t>(-2.f * KNOT_SCALE);
    z1[TFOIL] = vessl::cast<sample_t>(1.f * KNOT_SCALE);
    z2[TFOIL] = vessl::cast<sample_t>(0.f);

    static constexpr int LISSA = static_cast<int>(KnotType::LISSA);
    x1[LISSA] = vessl::cast<sample_t>(0.f);
    x2[LISSA] = vessl::cast<sample_t>(2.f * KNOT_SCALE);
    x3[LISSA] = vessl::phase_360; // TWO_PI;
    y1[LISSA] = vessl::cast<sample_t>(2.f * KNOT_SCALE);
    y2[LISSA] = 3 * vessl::cast<phasew_t>(vessl::phase_180); // 3*PI;
    y3[LISSA] = vessl::cast<sample_t>(0.f);
    z1[LISSA] = vessl::cast<sample_t>(0.f);
    z2[LISSA] = vessl::cast<sample_t>(1.f * KNOT_SCALE);
    
    // @todo TORUS scale is like 2x TFOIL and LISSA, try to fix that.
    // TORUS with c = 2 and a = 1:
    // x = (c + a*cos(p))*sin(q) = c * sin(q) + a * sin(q) * cos(p) => cx1 = 2, cx2 = sin(q), cx3 = 0
    // y = (c + a*cos(p))*cos(q) = c * cos(q) + a * cos(q) * cos(p) => cy1 = 2, cy2 = 0, cy3 = cos(q) 
    // z = a*sin(p) => cz1 = 0, cz2 = 1
    static constexpr int TORUS = static_cast<int>(KnotType::TORUS);
    x1[TORUS] = vessl::cast<sample_t>(2.f * KNOT_SCALE) * TORUS_SCALE;
    x2[TORUS] = vessl::cast<sample_t>(0.f); /*sin(qt)*/
    x3[TORUS] = vessl::phase_zero;
    y1[TORUS] = vessl::cast<sample_t>(2.f * KNOT_SCALE) * TORUS_SCALE;
    y2[TORUS] = vessl::phase_zero;
    y3[TORUS] = vessl::cast<sample_t>(0.f); /*cos(qt)*/
    z1[TORUS] = vessl::cast<sample_t>(0.f);
    // technically should be 1.f * KNOT_SCALE but that makes for a very flat torus
    z2[TORUS] = vessl::cast<sample_t>(6.f * KNOT_SCALE) * TORUS_SCALE;
  }

private:
  [[nodiscard]] VESSL_INLINE static coord_t sample
  (
    phase_t pt, phase_t qt, phase_t zt,
    sample_t cx1, sample_t cx2, phase_t cx3,
    sample_t cy1, phase_t cy2, sample_t cy3,
    sample_t cz1, sample_t cz2
  )
  {
    return coord_t(
      cx1 * vessl::math::sin<sample_t>(qt) + cx2 * vessl::math::cos<sample_t>(pt + cx3),
      cy1 * vessl::math::cos<sample_t>(qt + cy2) + cy3 * vessl::math::cos<sample_t>(pt),
      cz1 * vessl::math::sin<sample_t>(3 * zt) + cz2 * vessl::math::sin<sample_t>(pt)
    );
  }
  
public:
  [[nodiscard]] VESSL_INLINE parameter knotTypeA() const { return params.knotTypeA("knot type a", 'A'); }
  [[nodiscard]] VESSL_INLINE parameter knotTypeB() const { return params.knotTypeB("knot type b", 'B'); }
  // [0,1] sets morph amount from knot type a to knot type b
  [[nodiscard]] VESSL_INLINE parameter knotMorph() const { return params.knotMorph("knot morph", 'm'); }
  [[nodiscard]] VESSL_INLINE parameter knotP() const { return params.knotP("knot P", 'P'); }
  [[nodiscard]] VESSL_INLINE parameter knotQ() const { return params.knotQ("knot Q", 'Q'); }
  // frequency modulation of just the P part of the knot
  [[nodiscard]] VESSL_INLINE parameter knotModP() const { return params.knotModP("mod P amount", 'p'); }
  // frequency modulation of just the Q part of the knot
  [[nodiscard]] VESSL_INLINE parameter knotModQ() const { return params.knotModQ("mod Q amount", 'q'); }
  
  [[nodiscard]] VESSL_INLINE parameter frequency() const { return params.frequency("frequency", 'F'); }
  [[nodiscard]] VESSL_INLINE parameter phaseMod() const { return params.phaseMod("phase mod", 'f'); }

  [[nodiscard]] VESSL_INLINE const parameter_list& parameters() const override { return *this; }
  
  VESSL_INLINE coord_t generate() override
  {
    return generate<true>();
  }

  template<bool SmoothPQ = true>
  VESSL_INLINE coord_t generate()
  {
    // calculate coefficients based on knot type and morph settings
    if (knotA != params.knotTypeA.value 
     || knotB != params.knotTypeB.value 
     || knotM != params.knotMorph.value)
    {
      int i = static_cast<int>(params.knotTypeA.value);
      int j = static_cast<int>(params.knotTypeB.value);
      
      knotA = params.knotTypeA.value;
      knotB = params.knotTypeB.value;

      // rate-limit morph to prevent fizz when the value jumps.
      // note: there is a balance that has to be struck here for Bastl Citadel,
      // if kminc is too small, this can cause buffer underrun with large value swings.
      static constexpr phase_t kminc = vessl::phase_360 / (360*4); 
      knotM = params.knotMorph.value > knotM 
        ? (knotM + vessl::math::min(kminc, params.knotMorph.value - knotM))
        : (knotM - vessl::math::min(kminc, knotM - params.knotMorph.value));

      cx1 = vessl::math::lerp(x1[i], x1[j], knotM);
      cx3 = vessl::cast<phase_t>(vessl::math::lerp(x3[i], x3[j], knotM));
      cy1 = vessl::math::lerp(y1[i], y1[j], knotM);
      cy2 = vessl::cast<phase_t>(vessl::math::lerp(y2[i], y2[j], knotM));
      cz1 = vessl::math::lerp(z1[i], z1[j], knotM);
      cz2 = vessl::math::lerp(z2[i], z2[j], knotM);
    }

    int i = static_cast<int>(knotA);
    int j = static_cast<int>(knotB);
    phase_t m = knotM;

    phase_t fm = params.phaseMod.value;
    int32_t kp = static_cast<int32_t>(params.knotP.value);
    int32_t kq = static_cast<int32_t>(params.knotQ.value);

    // the four phases we need for sampling the curves
    // are calculated as multiples of phases running
    // at the same frequency as phaseZ (with phase modulation added).
    // this keeps the four curves properly aligned for blending.
    phase_t phaseP1 = phaseP * kp + fm;
    phase_t phaseQ1 = phaseQ * kq + fm;
    phase_t phaseZM = phaseZ + fm;
    phase_t phaseT1 = phaseQ1;
    
    x2[static_cast<int>(KnotType::TORUS)] = vessl::math::sin<sample_t>(phaseT1) * TORUS_SCALE;
    y3[static_cast<int>(KnotType::TORUS)] = vessl::math::cos<sample_t>(phaseT1) * TORUS_SCALE;

    sample_t cx2 = vessl::math::lerp(x2[i], x2[j], m);
    sample_t cy3 = vessl::math::lerp(y3[i], y3[j], m);

    coord_t a = sample(phaseP1, phaseQ1, phaseZM , cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

    // support fractional P and Q values by generating a curve
    // that is a bilinear interpolation of phase-sync'd curves
    // for F(P,Q), F(P+1,Q), F(P,Q+1), F(P+1,Q+1).
    if (SmoothPQ)
    {
      sample_t pd = vessl::cast<sample_t>(params.knotP.value - kp);
      sample_t qd = vessl::cast<sample_t>(params.knotQ.value - kq);
      
      //if (pd != 0 || qd != 0)
      {
        phase_t phaseP2 = phaseP1 + phaseP;
        phase_t phaseQ2 = phaseQ1 + phaseQ;
        phase_t phaseT2 = phaseQ2;

        coord_t b = sample(phaseP2, phaseQ1, phaseZM, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

        x2[static_cast<int>(KnotType::TORUS)] = vessl::math::sin<sample_t>(phaseT2) * TORUS_SCALE;
        y3[static_cast<int>(KnotType::TORUS)] = vessl::math::cos<sample_t>(phaseT2) * TORUS_SCALE;

        cx2 = vessl::math::lerp(x2[i], x2[j], m); // interp(x2, i, j, lerp);
        cy3 = vessl::math::lerp(y3[i], y3[j], m); // interp(y3, i, j, lerp);

        coord_t c = sample(phaseP1, phaseQ2, phaseZM, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);
        coord_t d = sample(phaseP2, phaseQ2, phaseZM, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

        // @todo return to pretty arithmetic statements when the operators are working correctly.
        auto aa = a.as_array();
        auto bb = b.as_array();
        auto cc = c.as_array();
        auto dd = d.as_array();
        //a = a + (b - a) * pd;
        coord_t a_b_p(b);
        a_b_p.as_array().subtract(aa).scale(pd).add(aa).copy_to(aa);
        //b = c + (d - c) * pd;
        coord_t c_d_p(d);
        c_d_p.as_array().subtract(cc).scale(pd).add(cc).copy_to(bb);
        //a = a + (b - a) * qd;
        bb.subtract(aa).scale(qd).add(aa).copy_to(aa);
      }
    }

    phase_t freqZ = dt * params.frequency.value;
    phase_t freqP = freqZ * (1.f+params.knotModP.value);
    phase_t freqQ = freqZ * (1.f+params.knotModQ.value);
    phaseP += freqP;
    phaseQ += freqQ;
    phaseZ += freqZ;

    xyz_ = a;
    return a;
  }
  
  // to help with debugging
  VESSL_INLINE phase_t pz() const { return phaseZ; }
  VESSL_INLINE phase_t pp() const { return phaseP; }
  VESSL_INLINE phase_t pq() const { return phaseQ; }
  VESSL_INLINE phase_t pi() const { return dt; }
  VESSL_INLINE float   sr() const { return sr_; }
  VESSL_INLINE coord_t xyz() const { return xyz_; }
  
protected:
  [[nodiscard]] VESSL_INLINE parameter element_at(vessl::size_t index) const override
  {
    parameter p[num] = { knotTypeA(), knotTypeB(), knotMorph(), knotP(), knotQ(), knotModP(), knotModQ(), frequency(), phaseMod() };
    return p[index];
  }

public:
  static KnotOscillator* create(float sr)
  {
    return new KnotOscillator(sr);
  }

  static void destroy(const KnotOscillator* knoscil)
  {
    delete knoscil;
  }
};
