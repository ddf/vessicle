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
  
  static constexpr int KNOT_TYPE_COUNT = static_cast<int>(KnotType::COUNT);
  
  using coord_t   = vessl::vector3<T>;
  using phase_t   = vessl::phase_t;
  using analog_t  = vessl::analog_t;
  using param     = vessl::parameter;
  using desc      = param::desc;
  using digital_p = vessl::digital_p;
  using analog_p  = vessl::analog_p;
  using knot_p    = vessl::param<KnotType>;
  using phase_p   = vessl::phase_p;
  
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
  analog_t x1[KNOT_TYPE_COUNT];
  analog_t x2[KNOT_TYPE_COUNT];
  uint64_t x3[KNOT_TYPE_COUNT];
  analog_t y1[KNOT_TYPE_COUNT];
  uint64_t y2[KNOT_TYPE_COUNT];
  analog_t y3[KNOT_TYPE_COUNT];
  analog_t z1[KNOT_TYPE_COUNT];
  analog_t z2[KNOT_TYPE_COUNT];
  
  phase_t phaseP;
  phase_t phaseQ;
  phase_t phaseZ;
  float   sampRate;

public:
  explicit KnotOscillator(float sampleRate)
    : params()
    , phaseP(vessl::PHASE_ZERO), phaseQ(vessl::PHASE_ZERO), phaseZ(vessl::PHASE_ZERO)
    , sampRate(sampleRate)
  {
    params.knotP.value = 1;
    params.knotQ.value = 1;
    params.knotTypeA.value = KnotType::TFOIL;
    params.knotTypeB.value = KnotType::LISSA;
    params.frequency.value = 1;
    
    // @todo use int for coefficients store radian coefficients as phase_t
    static constexpr int TFOIL = static_cast<int>(KnotType::TFOIL);
    x1[TFOIL] = 1;
    x2[TFOIL] = 2;
    x3[TFOIL] = 3*vessl::PHASE_HALF/2; // 3*PI/2;
    y1[TFOIL] = 1;
    y2[TFOIL] = 0;
    y3[TFOIL] = -2;
    z1[TFOIL] = 1;
    z2[TFOIL] = 0;

    static constexpr int LISSA = static_cast<int>(KnotType::LISSA);
    x1[LISSA] = 0;
    x2[LISSA] = 2;
    x3[LISSA] = vessl::PHASE_MAX; // TWO_PI;
    y1[LISSA] = 2;
    y2[LISSA] = 3*vessl::PHASE_HALF; // 3*PI;
    y3[LISSA] = 0;
    z1[LISSA] = 0;
    z2[LISSA] = 1;
    
    // TORUS with c = 2 and a = 1:
    // x = (c + a*cos(p))*sin(q) = c * sin(q) + a * sin(q) * cos(p) => cx1 = 2, cx2 = sin(q), cx3 = 0
    // y = (c + a*cos(p))*cos(q) = c * cos(q) + a * cos(q) * cos(p) => cy1 = 2, cy2 = 0, cy3 = cos(q) 
    // z = a*sin(p) => cz1 = 0, cz2 = 1
    static constexpr int TORUS = static_cast<int>(KnotType::TORUS);
    x1[TORUS] = 2;
    x2[TORUS] = 0; /*sin(qt)*/
    x3[TORUS] = 0;
    y1[TORUS] = 2;
    y2[TORUS] = 0;
    y3[TORUS] = 0; /*cos(qt)*/
    z1[TORUS] = 0;
    z2[TORUS] = 1;
  }

private:
  static coord_t sample(phase_t pt, phase_t qt, phase_t zt,
    analog_t cx1, analog_t cx2, phase_t cx3,
    analog_t cy1, phase_t  cy2, analog_t cy3,
    analog_t cz1, analog_t cz2)
  {
    return coord_t(
      cx1 * vessl::math::sinz<analog_t>(qt) + cx2 * vessl::math::cosz<analog_t>(pt + cx3),
      cy1 * vessl::math::cosz<analog_t>(qt + cy2) + cy3 * vessl::math::cosz<analog_t>(pt),
      cz1 * vessl::math::sinz<analog_t>(3 * zt) + cz2 * vessl::math::sinz<analog_t>(pt)
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
    
    phase_t  m = params.knotMorph.value;

    analog_t cx1 = vessl::easing::lerpp(x1[i], x1[j], m);
    phase_t  cx3 = vessl::cast<phase_t>(vessl::easing::lerpp(x3[i], x3[j], m));
    analog_t cy1 = vessl::easing::lerpp(y1[i], y1[j], m);
    phase_t  cy2 = vessl::cast<phase_t>(vessl::easing::lerpp(y2[i], y2[j], m));
    analog_t cz1 = vessl::easing::lerpp(z1[i], z1[j], m);
    analog_t cz2 = vessl::easing::lerpp(z2[i], z2[j], m);

    phase_t fm = params.phaseMod.value;
    analog_t kp = vessl::math::floor(params.knotP.value);
    analog_t kq = vessl::math::floor(params.knotQ.value);

    // the four phases we need for sampling the curves
    // are calculated as multiples of phases running
    // at the same frequency as phaseZ (with phase modulation added).
    // this keeps the four curves properly aligned for blending.
    phase_t phaseP1 = phaseP * static_cast<phase_t>(kp) + fm;
    phase_t phaseQ1 = phaseQ * static_cast<phase_t>(kq) + fm;
    phase_t phaseT1 = phaseQ1;
    
    x2[static_cast<int>(KnotType::TORUS)] = vessl::math::sinz<analog_t>(phaseT1);
    y3[static_cast<int>(KnotType::TORUS)] = vessl::math::cosz<analog_t>(phaseT1);

    T cx2 = vessl::easing::lerpp(x2[i], x2[j], m); // interp(x2, i, j, lerp);
    T cy3 = vessl::easing::lerpp(y3[i], y3[j], m); // interp(y3, i, j, lerp);

    coord_t a = sample(phaseP1, phaseQ1, phaseZ + fm, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

    // support fractional P and Q values by generating a curve
    // that is a bilinear interpolation of phase-sync'd curves
    // for F(P,Q), F(P+1,Q), F(P,Q+1), F(P+1,Q+1).
    if (smooth_pq)
    {
      analog_t pd = params.knotP.value - kp;
      analog_t qd = params.knotQ.value - kq;
      phase_t phaseP2 = phaseP * (static_cast<phase_t>(kp) + 1) + fm;
      phase_t phaseQ2 = phaseQ * (static_cast<phase_t>(kq) + 1) + fm;
      phase_t phaseT2 = phaseQ2;

      coord_t b = sample(phaseP2, phaseQ1, phaseZ + fm, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

      x2[static_cast<int>(KnotType::TORUS)] = vessl::math::sinz<analog_t>(phaseT2);
      y3[static_cast<int>(KnotType::TORUS)] = vessl::math::cosz<analog_t>(phaseT2);

      cx2 = vessl::easing::lerpp(x2[i], x2[j], m); // interp(x2, i, j, lerp);
      cy3 = vessl::easing::lerpp(y3[i], y3[j], m); // interp(y3, i, j, lerp);

      coord_t c = sample(phaseP1, phaseQ2, phaseZ + fm, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);
      coord_t d = sample(phaseP2, phaseQ2, phaseZ + fm, cx1, cx2, cx3, cy1, cy2, cy3, cz1, cz2);

      a = a + (b - a) * pd;
      b = c + (d - c) * pd;
      a = a + (b - a) * qd;
    }

    analog_t freqZ = params.frequency.value / sampRate;
    analog_t freqP = freqZ*(1+params.knotModP.value);
    analog_t freqQ = freqZ*(1+params.knotModQ.value);
    phaseP += vessl::cast<phase_t>(freqP);
    phaseQ += vessl::cast<phase_t>(freqQ);
    phaseZ += vessl::cast<phase_t>(freqZ);

    return a;
  }
  
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
