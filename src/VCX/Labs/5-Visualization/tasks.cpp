#include "Labs/5-Visualization/tasks.h"

#include <numbers>
#include <math.h>
#include <stdlib.h>

using VCX::Labs::Common::ImageRGB;
namespace VCX::Labs::Visualization {

    struct CoordinateStates {
        const float inf = 1e9;
        std::vector<Car> data;
        std::vector<bool> is_vaild;
        CoordinateStates(const std::vector<Car> &cars)
        {
            this -> data = cars;
            is_vaild.resize((int)(cars.size()));
            for(int i = 0; i < (int)(is_vaild.size()); ++i) is_vaild[i] = 1;
        }
        bool Update(InteractProxy const &proxy)
        {
            return true;
        }
        float get_value(Car &car, int id)
        {
            switch(id)
            {
                case 0:
                    return car.cylinders;
                case 1:
                    return car.displacement;
                case 2:
                    return car.weight;
                case 3:
                    return car.horsepower;
                case 4:
                    return car.acceleration;  
                case 5:
                    return car.mileage;
                case 6:
                    return car.year;
                default:
                    assert(0);
                    return 0;
            }
        }
        void Paint(Common::ImageRGB &input)
        {
            const int width = input.GetSizeX();
            const int height = input.GetSizeY();
            SetBackGround(input, glm::vec4(1.0f, 1.0f, 1.0f, 0.0f));
            std::vector<std::string> names = {
                "cylinders",
                "displacement",
                "weight",
                "horsepower",
                "acceleration",
                "mileage",
                "year"
            };
            ///calculate the maxs and mins
            std::vector<float> minv(names.size(), inf);
            std::vector<float> maxv(names.size(), -inf);
            for (auto &car : data) 
            {
                for(int i = 0; i < (int)(names.size()); ++i)
                {
                    minv[i] = std::min(minv[i], get_value(car, i));
                    maxv[i] = std::max(maxv[i], get_value(car, i));
                }
            }
            const float left_padding = 0.10f;
            const float right_padding = 0.10f;
            const float up_padding = 0.20f;
            const float down_padding = 0.05f;
            const float rect_width = 0.01f;
            const float line_distance = (1.0f - left_padding - right_padding) / (names.size() - 1);
            const int dirctor_nums = 10;
            const float dirctor_distance = (1.0f - up_padding - down_padding) / dirctor_nums;
            const float line_width = 2;
            const glm::vec4 text_color = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
            for (size_t i = 0; i < names.size(); ++i) {
                float from_x = left_padding + i * line_distance;/// left-up-x
                float from_y = up_padding; 
                float to_x = left_padding + i * line_distance;
                float to_y = 1.0f - down_padding;
                float middle_x = left_padding + i * line_distance;
                DrawFilledRect(input, 
                        glm::vec4(200 / 256.0f, 198 / 256.0f, 198 / 256.0f, 1.0f), 
                        glm::vec2(from_x - rect_width / 2, from_y),
                        glm::vec2(rect_width, 1.0f - down_padding - up_padding)
                    );
                DrawLine(input, 
                        glm::vec4(0.0f, 0.0f, 0.0f, 0.5f),
                        glm::vec2(from_x, from_y),
                        glm::vec2(to_x, to_y),
                        line_width
                    );
                PrintText(input, text_color, glm::vec2(middle_x, from_y - 0.08f), 0.02f, names[i]);
                PrintText(input, text_color, glm::vec2(middle_x, from_y - 0.02f), 0.015f, std::to_string((int)ceil(maxv[i])));
                PrintText(input, text_color, glm::vec2(middle_x, to_y + 0.02f), 0.015f, std::to_string((int)floor(minv[i])));
            }
            for (auto &car : data) 
            {
                glm::vec2 lastPoint;
                for (int i = 0; i < (int)names.size(); ++i) 
                {
                    if(!is_vaild[i]) continue;
                    float x = left_padding + i * line_distance;
                    float value = get_value(car, i);
                    float y = up_padding + (1.0f - up_padding - down_padding) * (value - minv[i]) / (maxv[i] - minv[i]);
                    glm::vec2 point(x, y);
                    if (i > 0) DrawLine(input, glm::vec4(0.0f, 0.5f, 0.8f, 1.0f), lastPoint, point, 1.0f);
                    lastPoint = point;
                }
            }
        }
        // your code here
    };
    bool PaintParallelCoordinates(Common::ImageRGB &input, InteractProxy const &proxy, std::vector<Car> const &data, bool force) {
        static CoordinateStates states(data); // initialize
        bool change = states.Update(proxy); // update according to user input
        if (! force && ! change) return false; // determine to skip repainting
        states.Paint(input); // visualize
        return true;
    }

    void Get_border(VectorField2D const &field, glm::vec2 &pos) 
    {
        float x = std::clamp(pos.x, 0.0f, (float)(field.size.first - 1));
        float y = std::clamp(pos.y, 0.0f, (float)(field.size.second - 1));
        pos = glm::vec2(x, y);
    }

    // relize LIC algorithm
    void LIC(Common::ImageRGB &output, Common::ImageRGB const &noise, VectorField2D const &field, int const &step) {
    
        int width = noise.GetSizeX();
        int height = noise.GetSizeY();
        // using square normalization
        for (int x = 0; x < height; ++x) 
        {
            for (int y = 0; y < width; ++y) 
            {
                glm::vec2 pos(x + 0.5f, y + 0.5f); ///get the central position
                glm::vec3 current_color = glm::vec3(0.0f, 0.0f, 0.0f);
                float weight_sum = 0.0f;

                ///this point
                weight_sum += 1;
                current_color += noise.At(x, y);
                glm::vec2 npos = pos;
                for(int i = 1; i <= step; ++i) // forward
                {
                    glm::vec2 dir = field.At(npos.x, npos.y);
                    float dx = 0, dy = 0;
                    if(dir.x > 0) dx = (floor(npos.x) + 1 - npos.x) / dir.x;
                    else if(dir.x < 0) dx = (npos.x - (ceil(npos.x) - 1)) / -dir.x;
                    if(dir.y > 0) dy = (floor(npos.y) + 1 - npos.y) / dir.y;
                    else if(dir.y < 0) dy = (npos.y - (ceil(npos.y) - 1)) / -dir.y;
                    float dt = std::min(dx, dy);
                    npos += dt * dir;
                    Get_border(field, npos);
                    float weight = 1.0f * (step - i) * (step - i) / step / step;
                    weight_sum += weight;
                    current_color += noise.At((int)npos.x, (int)npos.y) * weight;
                    //printf("(%d, %d)step%d : dt = %.2lf, dir = (%.2lf, %.2lf)npos = (%.2lf %.2lf)\n", \
                    //    x, y, i, dt, dir.x, dir.y, npos.x, npos.y);
                }
                npos = pos;
                for(int i = 1; i <= step; ++i) // forward
                {
                    glm::vec2 dir = -1.0f * field.At(npos.x, npos.y);
                    float dx = 0, dy = 0;
                    if(dir.x > 0) dx = (floor(npos.x) + 1 - npos.x) / dir.x;
                    else if(dir.x < 0) dx = (npos.x - (ceil(npos.x) - 1)) / -dir.x;
                    if(dir.y > 0) dy = (floor(npos.y) + 1 - npos.y) / dir.y;
                    else if(dir.y < 0) dy = (npos.y - (ceil(npos.y) - 1)) / -dir.y;
                    float dt = std::min(dx, dy);
                    npos += dt * dir;
                    Get_border(field, npos);
                    float weight = 1.0f * (step - i) * (step - i) / step / step;
                    weight_sum += weight;
                    current_color += noise.At((int)npos.x, (int)npos.y) * weight;
                }
                // normalize
                current_color /= weight_sum;
                output.At(x, y) = current_color; //result
            }
        }
    }

}; // namespace VCX::Labs::Visualization