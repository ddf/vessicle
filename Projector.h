#pragma once

#include "vessl/vessl.h"

template<typename T>
class Projector : public vessl::unitProcessor<vessl::vector3<T>, vessl::vector2<T>>
{

};