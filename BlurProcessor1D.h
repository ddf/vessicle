#pragma once

#include "vessl/vessl.h"
#include "CircularTexture.h"
#include "BlurKernel.h"

enum class BlurAxis : uint8_t
{
  X,
  Y
};

enum class TextureSizeType : uint8_t
{
  Integral,
  Fractional,
};

// note: TextureSizeType is used to specialize on float to support textures with non-integral dimensions
// probably there is a less confusing way to do this.
template<BlurAxis Axis, TextureSizeType TextureSizeType = TextureSizeType::Integral>
class BlurProcessor1D : public vessl::unit_processor<float>, protected vessl::plist<1>
{  
  using parameter = vessl::parameter;
  
  using size_t = vessl::size_t;
  using TextureType = CircularTexture<float>;
public:
  explicit BlurProcessor1D() = default;

  BlurProcessor1D(float sampleRate, float* textureData, size_t textureSizeX, size_t textureSizeY, const BlurKernel& kernel)
    : unit_processor(), texture(textureData, textureSizeX * textureSizeY, textureSizeX, textureSizeY)
    , kernel(kernel) 
  { params.textureSize.value = static_cast<float>(textureSizeX); }
  
  const parameter_list& parameters() const override { return *this; }

  parameter textureSize() const { return params.textureSize("Texture Size", 't'); }
  
  float process(const float& input) override
  {
    texture.write(input);

    if (TextureSizeType == TextureSizeType::Integral)
    {
      return processIntegral();
    }

    return processFractional();
  }

  using unit_processor::process;
  
  static BlurProcessor1D* create(float sampleRate, size_t maxTextureSize, const BlurKernel& blurKernel)
  {
    // HACK: extra memory to ensure we don't read outside the bounds
    maxTextureSize += 2;
    if (Axis == BlurAxis::X)
    {
      return new BlurProcessor1D(sampleRate, new float[maxTextureSize], maxTextureSize, 1, blurKernel);
    }

    return new BlurProcessor1D(sampleRate, new float[maxTextureSize*maxTextureSize], maxTextureSize, maxTextureSize, blurKernel);
  }

  static void destroy(BlurProcessor1D* blur)
  {
    delete[] blur->texture.getData();
    delete blur;
  }
  
protected:
  parameter element_at(vessl::size_t index) const override
  {
    parameter p[num] = { textureSize() }; return p[index];
  }

private:
  struct
  {
    vessl::analog_p textureSize;
  } params;
  TextureType texture;
  BlurKernel kernel;
  
  float processIntegral()
  {
    float v = 0;
    float c = kernel.getBlurSize() * 0.5f;
    size_t texSize = textureSize().read_digital();
    if (Axis == BlurAxis::X)
    {
      TextureType tex = texture.subtexture(texSize, 1);
      for (const BlurKernelSample& samp : kernel)
      {
        v += tex.readBilinear(c + samp.offset, 0) * samp.weight;
      }
    }
    else
    {
      TextureType tex = texture.subtexture(texSize, texSize);
      for (const BlurKernelSample& samp : kernel)
      {
        v += tex.readBilinear(0, c + samp.offset) * samp.weight;
      }
    }

    return v;
  }

  float processFractional()
  {
    float v = 0;
    float c = kernel.getBlurSize() * 0.5f;
    float texSize = textureSize().read_analog();
    if (Axis == BlurAxis::X)
    {
      for (const BlurKernelSample& samp : kernel)
      {
        //v += texture.readBilinear(c + samp.offset, 0) * samp.weight;

        float x = (c + samp.offset) * texSize;
        int x1 = static_cast<int>(x);
        int x2 = x1 + 1;
        float xt = x - static_cast<float>(x1);

        v += vessl::math::lerp(texture.read(x1, 0), texture.read(x2, 0), xt) * samp.weight;
      }
    }
    else
    {
      size_t texSizeLow = (size_t)texSize;
      size_t texSizeHi = texSizeLow + 1;
      float texSizeBlend = texSize - static_cast<float>(texSizeLow);
      // read from our two integral texture sizes centered on the same position in the big texture
      // adding texSize to the blur-based offset prevents reading past the write head,
      // which introduces a delay-like echo.
      float readOffset = (texSize) * (texSize) * c + texSize;
      int x1 = static_cast<int>(readOffset);
      int x2 = x1 + 1;
      float xt = readOffset - static_cast<float>(x1);
      TextureType texA = texture.subtexture(texSizeLow, texSizeLow);
      TextureType texB = texture.subtexture(texSizeHi, texSizeHi);
      for (const BlurKernelSample& samp : kernel)
      {
        //v += Interpolator::linear(texA.readBilinear(u1, coord), texB.readBilinear(u2, coord), texSizeBlend) * samp.weight;

        // this is essentially the same as the line commented out above,
        // but instead of generating two different u coordinates that align the reads,
        // we just calculate the x offset, which will be bigger than the texture dimension,
        // so the y coordinate can be +/- around that.
        float ya = samp.offset * static_cast<float>(texSizeLow);
        int ya1 = ya < 0 ? static_cast<int>(ya) - 1 : static_cast<int>(ya);
        int ya2 = ya1 + 1;
        float yat = ya - static_cast<float>(ya1);

        float yb = samp.offset * static_cast<float>(texSizeHi);
        int yb1 = yb < 0 ? static_cast<int>(yb) - 1 : static_cast<int>(yb);
        int yb2 = yb1 + 1;
        float ybt = yb - static_cast<float>(yb1);

        float xa1 = vessl::math::lerp(texA.read(x1, ya1), texA.read(x2, ya1), xt);
        float xa2 = vessl::math::lerp(texA.read(x1, ya2), texA.read(x2, ya2), xt);
        float va  = vessl::math::lerp(xa1, xa2, yat);

        float xb1 = vessl::math::lerp(texB.read(x1, yb1), texB.read(x2, yb1), xt);
        float xb2 = vessl::math::lerp(texB.read(x1, yb2), texB.read(x2, yb2), xt);
        float vb  = vessl::math::lerp(xb1, xb2, ybt);

        v += vessl::math::lerp(va, vb, texSizeBlend) * samp.weight;
      }
    }

    return v;
  }
};