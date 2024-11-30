#include "Labs/3-Rendering/tasks.h"

namespace VCX::Labs::Rendering {

    glm::vec4 GetTexture(Engine::Texture2D<Engine::Formats::RGBA8> const & texture, glm::vec2 const & uvCoord) {
        if (texture.GetSizeX() == 1 || texture.GetSizeY() == 1) return texture.At(0, 0);
        glm::vec2 uv      = glm::fract(uvCoord);
        uv.x              = uv.x * texture.GetSizeX() - .5f;
        uv.y              = uv.y * texture.GetSizeY() - .5f;
        std::size_t xmin  = std::size_t(glm::floor(uv.x) + texture.GetSizeX()) % texture.GetSizeX();
        std::size_t ymin  = std::size_t(glm::floor(uv.y) + texture.GetSizeY()) % texture.GetSizeY();
        std::size_t xmax  = (xmin + 1) % texture.GetSizeX();
        std::size_t ymax  = (ymin + 1) % texture.GetSizeY();
        float       xfrac = glm::fract(uv.x), yfrac = glm::fract(uv.y);
        return glm::mix(glm::mix(texture.At(xmin, ymin), texture.At(xmin, ymax), yfrac), glm::mix(texture.At(xmax, ymin), texture.At(xmax, ymax), yfrac), xfrac);
    }

    glm::vec4 GetAlbedo(Engine::Material const & material, glm::vec2 const & uvCoord) {
        glm::vec4 albedo       = GetTexture(material.Albedo, uvCoord);
        glm::vec3 diffuseColor = albedo;
        return glm::vec4(glm::pow(diffuseColor, glm::vec3(2.2)), albedo.w);
    }

    /******************* 1. Ray-triangle intersection *****************/
    bool IntersectTriangle(Intersection & output, Ray const & ray, glm::vec3 const & p1, glm::vec3 const & p2, glm::vec3 const & p3) {
        // your code here
        glm::vec3 D = ray.Direction;
        glm::vec3 E1 = p2 - p1, E2 = p3 - p1, T = ray.Origin - p1;
        float d = dot(cross(D, E2), E1);
        float t = dot(cross(T, E1), E2);
        float u = dot(cross(D, E2), T);
        float v = dot(cross(T, E1), D);
        if(fabs(d) < 1e-9) return false;
        t /= d, u /= d, v /= d;
        if(u >= 0 && v >= 0 && u + v <= 1) 
        {
            output = {t, u, v};
            return true;
        }
        return false;
    }
    //阅读并填补 tasks.cpp 中光线追踪函数 RayTrace 中的着色部分。
    //光线与场景求交形成的位置、法向、反照率、吸收率、透明度、高光衰减指数都已在函数内给出。在
    //完成着色部分之后，阅读讲义中 Shadow Ray 的思想，最后在光线追踪中实现阴影。你可以使用
    //auto hit = intersector.IntersectRay(Ray(pos, dir));
    //语句来对光线求交；为了简化实现，在计算阴影过程中，如果光源与着色点之间存在遮挡物，
    //近似认为 alpha < 0.2 的遮挡物视为透明，而 alpha >= 0.2 的遮挡物视为不可透过。
    //这里透明度 alpha = hit.IntersectAlbedo.w 。

    glm::vec3 RayTrace(const RayIntersector & intersector, Ray ray, int maxDepth, bool enableShadow) {
        glm::vec3 color(0.0f);
        glm::vec3 weight(1.0f);
        //1. 对于每一个像素点 (x, y)，从观察者的眼睛发出一条穿过 (x, y) 的光线．
        //2. 找出这条光线与场景中的物体第一次发生相交的位置．
        for (int depth = 0; depth < maxDepth; depth++) {
            auto rayHit = intersector.IntersectRay(ray);
            if (! rayHit.IntersectState) return color;
            const glm::vec3 pos       = rayHit.IntersectPosition;//交点
            const glm::vec3 n         = rayHit.IntersectNormal;//法向
            const glm::vec3 kd        = rayHit.IntersectAlbedo;//反照率
            const glm::vec3 ks        = rayHit.IntersectMetaSpec;//吸收率
            const float     alpha     = rayHit.IntersectAlbedo.w;//透明度
            const float     shininess = rayHit.IntersectMetaSpec.w * 256;//衰减指数

            glm::vec3 result(0.0f);
            /******************* 2. Whitted-style ray tracing *****************/
            // your code here

            for (const Engine::Light & light : intersector.InternalScene->Lights) {
                glm::vec3 l;
                float     attenuation;
                /******************* 3. Shadow ray *****************/
                if (light.Type == Engine::LightType::Point) {
                    l           = light.Position - pos;
                    attenuation = 1.0f / glm::dot(l, l);//衰减(点光源)
                    if (enableShadow) {
                        auto shadow_hit = intersector.IntersectRay(Ray(pos, l));
                        if(shadow_hit.IntersectState && shadow_hit.IntersectAlbedo.w >= 0.2)
                        {
                            glm::vec3 temp = shadow_hit.IntersectPosition;
                            bool flag1 = temp.x <= pos.x;
                            bool flag2 = temp.x <= light.Position.x;
                            if((flag1 && !flag2) || (!flag1 && flag2)) continue;
                        }
                        // your code here
                    }
                } else if (light.Type == Engine::LightType::Directional) {
                    l           = light.Direction;
                    attenuation = 1.0f;
                    if (enableShadow) {
                        // your code here
                        auto shadow_hit = intersector.IntersectRay(Ray(pos, l));
                        if(shadow_hit.IntersectState && shadow_hit.IntersectAlbedo.w >= 0.2)
                        {
                            continue ;
                        }
                    }
                }

                /******************* 2. Whitted-style ray tracing *****************/
                // your code here
                l = normalize(l);
                glm::vec3 ia = intersector.InternalScene->AmbientIntensity;
                glm::vec3 id = attenuation * light.Intensity;
                glm::vec3 h = normalize((l + ray.Direction));
                float spectacular = glm::pow(glm::max(glm::dot(n, h), .0f),shininess);//高光
                float diffuse = glm::max(glm::dot(n, l), .0f);
                result += kd * (ia + id * diffuse) + ks * id * spectacular;
            }

            if (alpha < 0.9) {
                // refraction
                // accumulate color
                glm::vec3 R = alpha * glm::vec3(1.0f);
                color += weight * R * result;
                weight *= glm::vec3(1.0f) - R;

                // generate new ray
                ray = Ray(pos, ray.Direction);
            } else {
                // reflection
                // accumulate color
                glm::vec3 R = ks * glm::vec3(0.5f);
                color += weight * (glm::vec3(1.0f) - R) * result;
                weight *= R;

                // generate new ray
                glm::vec3 out_dir = ray.Direction - glm::vec3(2.0f) * n * glm::dot(n, ray.Direction);
                ray               = Ray(pos, out_dir);
            }
        }

        return color;
    }
} // namespace VCX::Labs::Rendering