#pragma once

#include "vessl/vessl.h"

// projects a 3D coordinate to 2D using the simplest version of the formula.
// @todo add additionl parameters like FOV?
template<typename T>
class Projector : public vessl::unitProcessor<vessl::frame::channels<T, 3>, vessl::frame::channels<T, 2>>
  , protected vessl::plist<1>
{
  using analog_t = vessl::analog_t;
  using analog_p = vessl::analog_p;

  struct
  {
    analog_p zoom;
  } params;
  
public:
  using input_t = vessl::frame::channels<T, 3>;
  using output_t = vessl::frame::channels<T, 2>;
  using param = vessl::parameter;
  
  explicit Projector(analog_t initialZoom = 1.0f)
  {
    params.zoom.value = initialZoom;
  }
  
  [[nodiscard]] VESSL_INLINE param zoom() const{ return params.zoom({"zoom", 'Z', analog_p::type}); }
  
  const parameters& getParameters() const override { return *this; }
  
  VESSL_INLINE output_t process(const input_t& input) override
  {
    analog_t zm = params.zoom.value;
    analog_t cz = vessl::cast<analog_t>(input.z());
    T proj = vessl::cast<T>(1.f / (cz+zm));
    return output_t(input.x()*proj, input.y()*proj);
  }
  
  static Projector* create(analog_t initialZoom = 1.0f)
  {
    return new Projector(initialZoom);
  }
  
  static void destroy(const Projector* projector)
  {
    delete projector;
  }

protected:
  vessl::parameter elementAt(vessl::size_t index) const override
  {
    VASSERT(index == 0, "Attempted to access invalid index in Projector::elementAt");
    return zoom();
  }
};