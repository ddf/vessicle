#pragma once

#include "vessl/vessl.h"

template<typename T = vessl::analog_t>
class KnotOscillator : public vessl::unitGenerator<vessl::vector3<T>>
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

  using sample_t  = T;
  using coord_t   = vessl::vector3<sample_t>;
  using phase_t   = vessl::phase_t;
  using phasew_t  = vessl::size_t;
  using analog_t  = vessl::analog_t;
  using param     = vessl::parameter;
  using desc      = param::desc;
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
  coord_t a;

public:
  explicit KnotOscillator(float sampleRate)
    : params()
    , phaseP(vessl::PHASE_ZERO), phaseQ(vessl::PHASE_ZERO), phaseZ(vessl::PHASE_ZERO)
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
    x3[TFOIL] = 3 * vessl::cast<phasew_t>(vessl::PHASE_90); // 3*PI/2;
    y1[TFOIL] = vessl::cast<sample_t>(1.f * KNOT_SCALE);
    y2[TFOIL] = vessl::PHASE_ZERO;
    y3[TFOIL] = vessl::cast<sample_t>(-2.f * KNOT_SCALE);
    z1[TFOIL] = vessl::cast<sample_t>(1.f * KNOT_SCALE);
    z2[TFOIL] = vessl::cast<sample_t>(0.f);

    static constexpr int LISSA = static_cast<int>(KnotType::LISSA);
    x1[LISSA] = vessl::cast<sample_t>(0.f);
    x2[LISSA] = vessl::cast<sample_t>(2.f * KNOT_SCALE);
    x3[LISSA] = vessl::PHASE_360; // TWO_PI;
    y1[LISSA] = vessl::cast<sample_t>(2.f * KNOT_SCALE);
    y2[LISSA] = 3 * vessl::cast<phasew_t>(vessl::PHASE_180); // 3*PI;
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
    x3[TORUS] = vessl::PHASE_ZERO;
    y1[TORUS] = vessl::cast<sample_t>(2.f * KNOT_SCALE) * TORUS_SCALE;
    y2[TORUS] = vessl::PHASE_ZERO;
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
  [[nodiscard]] VESSL_INLINE param knotTypeA() const { return params.knotTypeA({ "knot type a", 'A', knot_p::type }); }
  [[nodiscard]] VESSL_INLINE param knotTypeB() const { return params.knotTypeB({ "knot type b", 'B', knot_p::type }); }
  // [0,1] sets morph amount from knot type a to knot type b
  [[nodiscard]] VESSL_INLINE param knotMorph() const { return params.knotMorph({ "knot morph", 'm', phase_p::type }); }
  [[nodiscard]] VESSL_INLINE param knotP() const { return params.knotP({ "knot P", 'P', analog_p::type }); }
  [[nodiscard]] VESSL_INLINE param knotQ() const { return params.knotQ({ "knot Q", 'Q', analog_p::type }); }
  // frequency modulation of just the P part of the knot
  [[nodiscard]] VESSL_INLINE param knotModP() const { return params.knotModP({ "mod P amount", 'p', analog_p::type }); }
  // frequency modulation of just the Q part of the knot
  [[nodiscard]] VESSL_INLINE param knotModQ() const { return params.knotModQ({ "mod Q amount", 'q', analog_p::type }); }
  
  [[nodiscard]] VESSL_INLINE param frequency() const { return params.frequency({ "frequency", 'F', analog_p::type }); }
  [[nodiscard]] VESSL_INLINE param phaseMod() const { return params.phaseMod({ "phase mod", 'f', phase_p::type }); }

  [[nodiscard]] VESSL_INLINE const vessl::parameters& getParameters() const override { return *this; }
  
  VESSL_INLINE coord_t generate() override
  {
    return generate<true>();
  }

  template<bool smooth_pq = true>
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
      knotM = params.knotMorph.value;

      cx1 = vessl::easing::lerpp(x1[i], x1[j], knotM);
      cx3 = vessl::cast<phase_t>(vessl::easing::lerpp(x3[i], x3[j], knotM));
      cy1 = vessl::easing::lerpp(y1[i], y1[j], knotM);
      cy2 = vessl::cast<phase_t>(vessl::easing::lerpp(y2[i], y2[j], knotM));
      cz1 = vessl::easing::lerpp(z1[i], z1[j], knotM);
      cz2 = vessl::easing::lerpp(z2[i], z2[j], knotM);
    }

    int i = static_cast<int>(knotA);
    int j = static_cast<int>(knotB);
    phase_t m = knotM;

    phase_t fm = params.phaseMod.value;
    int32_t kp = (int32_t)(params.knotP.value);
    int32_t kq = (int32_t)(params.knotQ.value);

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

    sample_t cx2 = vessl::easing::lerpp(x2[i], x2[j], m); // interp(x2, i, j, lerp);
    sample_t cy3 = vessl::easing::lerpp(y3[i], y3[j], m); // interp(y3, i, j, lerp);

    a = sample(phaseP1, phaseQ1, phaseZM , cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

    // support fractional P and Q values by generating a curve
    // that is a bilinear interpolation of phase-sync'd curves
    // for F(P,Q), F(P+1,Q), F(P,Q+1), F(P+1,Q+1).
    if (smooth_pq)
    {
      sample_t pd = vessl::cast<T>(params.knotP.value - kp);
      sample_t qd = vessl::cast<T>(params.knotQ.value - kq);
      phase_t phaseP2 = phaseP * (kp + 1) + fm;
      phase_t phaseQ2 = phaseQ * (kq + 1) + fm;
      phase_t phaseT2 = phaseQ2;

      coord_t b = sample(phaseP2, phaseQ1, phaseZM, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

      x2[static_cast<int>(KnotType::TORUS)] = vessl::math::sin<sample_t>(phaseT2) * TORUS_SCALE;
      y3[static_cast<int>(KnotType::TORUS)] = vessl::math::cos<sample_t>(phaseT2) * TORUS_SCALE;

      cx2 = vessl::easing::lerpp(x2[i], x2[j], m); // interp(x2, i, j, lerp);
      cy3 = vessl::easing::lerpp(y3[i], y3[j], m); // interp(y3, i, j, lerp);

      coord_t c = sample(phaseP1, phaseQ2, phaseZM, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);
      coord_t d = sample(phaseP2, phaseQ2, phaseZM, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

      a = a + (b - a) * pd;
      b = c + (d - c) * pd;
      a = a + (b - a) * qd;
    }

    phase_t freqZ = dt * params.frequency.value;
    phase_t freqP = freqZ; //(1.f+params.knotModP.value);
    phase_t freqQ = freqZ; //(1.f+params.knotModQ.value);
    phaseP += freqP;
    phaseQ += freqQ;
    phaseZ += freqZ;

    return a;
  }
  
  // to help with debugging
  VESSL_INLINE phase_t pz() const { return phaseZ; }
  VESSL_INLINE phase_t pp() const { return phaseP; }
  VESSL_INLINE phase_t pq() const { return phaseQ; }
  VESSL_INLINE phase_t pi() const { return dt; }
  VESSL_INLINE float   sr() const { return sr_; }
  VESSL_INLINE coord_t xyz() const { return a; }
  
protected:
  [[nodiscard]] VESSL_INLINE param elementAt(vessl::size_t index) const override
  {
    param p[plsz] = { knotTypeA(), knotTypeB(), knotMorph(), knotP(), knotQ(), knotModP(), knotModQ(), frequency(), phaseMod() };
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
