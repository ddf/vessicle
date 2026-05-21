#pragma once

#include "vessl/vessl.h"

// projects a 3D coordinate to 2D using the simplest version of the formula.
// @todo add additionl parameters like FOV?
template<typename T>
class Projector : public vessl::unit_processor<vessl::frame::channels<T, 3>, vessl::frame::channels<T, 2>>
  , protected vessl::plist<1>
{
public:
  using analog_t = vessl::analog_t;
  using analog_p = vessl::analog_p;
  using input_t = vessl::frame::channels<T, 3>;
  using output_t = vessl::frame::channels<T, 2>;
  using param = vessl::parameter;
  
  explicit Projector(analog_t initial_zoom = 1.0f)
  {
    params_.zoom.value = initial_zoom;
  }
  
  [[nodiscard]] VESSL_INLINE param zoom() const{ return params_.zoom({"zoom", 'Z', analog_p::type}); }
  
  const parameter_list& parameters() const override { return *this; }
  
  VESSL_INLINE output_t process(const input_t& input) override
  {
    analog_t zm = params_.zoom.value;
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
  vessl::parameter element_at(vessl::size_t index) const override
  {
    VASSERT(index == 0, "Attempted to access invalid index in Projector::elementAt");
    return zoom();
  }
  
private:
  struct
  {
    analog_p zoom;
  } params_;
};