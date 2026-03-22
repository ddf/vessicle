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

  using coord_t   = vessl::vector3<T>;
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
  T x1[KNOT_TYPE_COUNT];
  T x2[KNOT_TYPE_COUNT];
  phasew_t x3[KNOT_TYPE_COUNT];
  T y1[KNOT_TYPE_COUNT];
  phasew_t y2[KNOT_TYPE_COUNT];
  T y3[KNOT_TYPE_COUNT];
  T z1[KNOT_TYPE_COUNT];
  T z2[KNOT_TYPE_COUNT];
  
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
    x1[TFOIL] = 1.f * KNOT_SCALE;
    x2[TFOIL] = 2.f * KNOT_SCALE;
    x3[TFOIL] = 3 * vessl::cast<phasew_t>(vessl::PHASE_90); // 3*PI/2;
    y1[TFOIL] = 1.f * KNOT_SCALE;
    y2[TFOIL] = 0.f;
    y3[TFOIL] = -2.f * KNOT_SCALE;
    z1[TFOIL] = 1.f * KNOT_SCALE;
    z2[TFOIL] = 0.f;

    static constexpr int LISSA = static_cast<int>(KnotType::LISSA);
    x1[LISSA] = 0.f;
    x2[LISSA] = 2.f * KNOT_SCALE;
    x3[LISSA] = vessl::PHASE_360; // TWO_PI;
    y1[LISSA] = 2.f * KNOT_SCALE;
    y2[LISSA] = 3 * vessl::cast<phasew_t>(vessl::PHASE_180); // 3*PI;
    y3[LISSA] = 0.f;
    z1[LISSA] = 0.f;
    z2[LISSA] = 1.f * KNOT_SCALE;
    
    // @todo TORUS scale is like 2x TFOIL and LISSA, try to fix that.
    // TORUS with c = 2 and a = 1:
    // x = (c + a*cos(p))*sin(q) = c * sin(q) + a * sin(q) * cos(p) => cx1 = 2, cx2 = sin(q), cx3 = 0
    // y = (c + a*cos(p))*cos(q) = c * cos(q) + a * cos(q) * cos(p) => cy1 = 2, cy2 = 0, cy3 = cos(q) 
    // z = a*sin(p) => cz1 = 0, cz2 = 1
    static constexpr int TORUS = static_cast<int>(KnotType::TORUS);
    x1[TORUS] = 2.f * KNOT_SCALE;
    x2[TORUS] = 0.f; /*sin(qt)*/
    x3[TORUS] = 0.f;
    y1[TORUS] = 2.f * KNOT_SCALE;
    y2[TORUS] = 0.f;
    y3[TORUS] = 0.f; /*cos(qt)*/
    z1[TORUS] = 0.f;
    z2[TORUS] = 1.f * KNOT_SCALE;
  }

private:
  static coord_t sample(phase_t pt, phase_t qt, phase_t zt,
    T cx1, T cx2, phase_t cx3,
    T cy1, phase_t cy2, T cy3,
    T cz1, T cz2)
  {
    return coord_t(
      cx1 * vessl::math::sin<T>(qt) + cx2 * vessl::math::cos<T>(pt + cx3),
      cy1 * vessl::math::cos<T>(qt + cy2) + cy3 * vessl::math::cos<T>(pt),
      cz1 * vessl::math::sin<T>(zt + zt + zt) + cz2 * vessl::math::sin<T>(pt)
    );
  }
  
public:
  param knotTypeA() const { return params.knotTypeA({ "knot type a", 'A', knot_p::type }); }
  param knotTypeB() const { return params.knotTypeB({ "knot type b", 'B', knot_p::type }); }
  // [0,1] sets morph amount from knot type a to knot type b
  param knotMorph() const { return params.knotMorph({ "knot morph", 'm', phase_p::type }); }
  param knotP() const { return params.knotP({ "knot P", 'P', analog_p::type }); }
  param knotQ() const { return params.knotQ({ "knot Q", 'Q', analog_p::type }); }
  // frequency modulation of just the P part of the knot
  param knotModP() const { return params.knotModP({ "mod P amount", 'p', analog_p::type }); }
  // frequency modulation of just the Q part of the knot
  param knotModQ() const { return params.knotModQ({ "mod Q amount", 'q', analog_p::type }); }
  
  param frequency() const { return params.frequency({ "frequency", 'F', analog_p::type }); }
  param phaseMod() const { return params.phaseMod({ "phase mod", 'f', phase_p::type }); }

  [[nodiscard]] const vessl::parameters& getParameters() const override { return *this; }
  
  coord_t generate() override
  {
    return generate<true>();
  }

  template<bool smooth_pq = true>
  coord_t generate()
  {
    // calculate coefficients based on knot type and morph settings
    int i = static_cast<int>(params.knotTypeA.value);
    int j = static_cast<int>(params.knotTypeB.value);
    
    phase_t m = params.knotMorph.value;

    T cx1 = vessl::easing::lerpp(x1[i], x1[j], m);
    phase_t cx3 = vessl::cast<phase_t>(vessl::easing::lerpp(x3[i], x3[j], m));
    T cy1 = vessl::easing::lerpp(y1[i], y1[j], m);
    phase_t cy2 = vessl::cast<phase_t>(vessl::easing::lerpp(y2[i], y2[j], m));
    T cz1 = vessl::easing::lerpp(z1[i], z1[j], m);
    T cz2 = vessl::easing::lerpp(z2[i], z2[j], m);

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
    
    x2[static_cast<int>(KnotType::TORUS)] = vessl::math::sin<T>(phaseT1);
    y3[static_cast<int>(KnotType::TORUS)] = vessl::math::cos<T>(phaseT1);

    T cx2 = x2[i]; // vessl::easing::lerpp(x2[i], x2[j], m); // interp(x2, i, j, lerp);
    T cy3 = y3[i]; // vessl::easing::lerpp(y3[i], y3[j], m); // interp(y3, i, j, lerp);

    a = sample(phaseP1, phaseQ1, phaseZM , cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

    // support fractional P and Q values by generating a curve
    // that is a bilinear interpolation of phase-sync'd curves
    // for F(P,Q), F(P+1,Q), F(P,Q+1), F(P+1,Q+1).
    if (smooth_pq)
    {
      T pd = vessl::cast<T>(params.knotP.value - kp);
      T qd = vessl::cast<T>(params.knotQ.value - kq);
      phase_t phaseP2 = phaseP * (kp + 1) + fm;
      phase_t phaseQ2 = phaseQ * (kq + 1) + fm;
      phase_t phaseT2 = phaseQ2;

      coord_t b = sample(phaseP2, phaseQ1, phaseZM, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

      x2[static_cast<int>(KnotType::TORUS)] = vessl::math::sin<T>(phaseT2);
      y3[static_cast<int>(KnotType::TORUS)] = vessl::math::cos<T>(phaseT2);

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
  phase_t pz() const { return phaseZ; }
  phase_t pp() const { return phaseP; }
  phase_t pq() const { return phaseQ; }
  phase_t pi() const { return dt; }
  float   sr() const { return sr_; }
  coord_t xyz() const { return a; }
  
protected:
  [[nodiscard]] param elementAt(vessl::size_t index) const override
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
