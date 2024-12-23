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
        std::vector<std::string> names = {
                "cylinders",
                "displacement",
                "weight",
                "horsepower",
                "acceleration",
                "mileage",
                "year"
            };
        std::vector<float> minv;
        std::vector<float> maxv;
        std::vector<std::pair<float, float>> clicked;
        std::vector<int> type;
        const float left_padding = 0.10f;
        const float right_padding = 0.10f;
        const float up_padding = 0.20f;
        const float down_padding = 0.05f;
        const float rect_width = 0.02f;
        const float line_distance = (1.0f - left_padding - right_padding) / (names.size() - 1);
        const int dirctor_nums = 10;
        const float dirctor_distance = (1.0f - up_padding - down_padding) / dirctor_nums;
        const float line_width = 2;
        const glm::vec4 text_color = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
        CoordinateStates(const std::vector<Car> &cars)
        {
            this -> data = cars;
            is_vaild.resize((int)(cars.size()));
            for(int i = 0; i < (int)(is_vaild.size()); ++i) is_vaild[i] = 1;
            minv.resize(names.size());
            maxv.resize(names.size());
            for(int i = 0; i < (int)(names.size()); ++i) minv[i] = inf, maxv[i] = -inf;
            for(int i = 0; i < (int)(names.size()); ++i) clicked.push_back({0.0f, 1.0f});
            for(int i = 0; i < (int)(names.size()); ++i) type.push_back(0);
        }
        int find_block(float &x, float &y)
        {
            int block = -1;
            for(int i = 0; i < (int)(names.size()); ++i)
            {
                /// a more flexible width
                float from_x = left_padding + i * line_distance - line_distance / 2 ;/// left-up-x
                float from_y = up_padding; 
                float to_x = left_padding + i * line_distance + line_distance / 2;
                float to_y = 1.0f - down_padding;
                //printf("%.2f %.2f\n", x, y);
                if(x >= from_x && x <= to_x && y >= from_y && y <= to_y)
                {
                    block = i;
                    //printf("%d\n", i);
                    break;
                }
            }
            return block;
        }
        bool Update(InteractProxy const &proxy)
        {
            int block = -1;
            if(proxy.IsHovering()) 
            {
                float x = proxy.MousePos().x;
                float y = proxy.MousePos().y;
                block = find_block(x, y);
                if(block != -1) ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
            }
            if(proxy.IsClicking() && proxy.IsHovering()) 
            {
                float x = proxy.MousePos().x;
                float y = proxy.MousePos().y;
                if (y >= 0.03f && y <= 0.09f && x >= 0.7f && x <= 0.9f) 
                {
                    //puts("1");
                    for(int i = 0; i < (int)(is_vaild.size()); ++i) is_vaild[i] = 1;
                    for(int i = 0; i < (int)(names.size()); ++i) type[i] = 0;
                    for(int i = 0; i < (int)(names.size()); ++i) clicked[i] = {0.0f, 1.0f};
                }
                if(block != -1)
                {
                    clicked[block] = {y - 0.002f, y + 0.002f};
                    type[block] = 1;
                    for(int i = 0; i < (int)(is_vaild.size()); ++i) is_vaild[i] = 1;
                    for(int i = 0; i < (int)data.size(); ++i)
                    {
                        for(int j = 0; j < (int)(names.size()); ++j)
                        {
                            float value = get_value(data[i], j);
                            float yy = up_padding + (1.0f - up_padding - down_padding) * (value - minv[j]) / (maxv[j] - minv[j]);
                            if(!(clicked[j].first <= yy && yy <= clicked[j].second)) is_vaild[i] = 0;
                        }
                    }
                }
                return true;
            }
            if(proxy.IsDragging() && proxy.IsHovering() && block != -1)
            {
                float x = proxy.MousePos().x;
                float y = proxy.MousePos().y;
                float old_x = proxy.DraggingStartPoint().x;
                float old_y = proxy.DraggingStartPoint().y;
                int old_block = find_block(old_x, old_y);
                if(old_block == block) 
                {
                    clicked[block] = {std::min(y, old_y), std::max(y, old_y)};
                    type[block] = 2;
                    for(int i = 0; i < (int)(is_vaild.size()); ++i) is_vaild[i] = 1;
                    for(int i = 0; i < (int)data.size(); ++i)
                    {
                        for(int j = 0; j < (int)(names.size()); ++j)
                        {
                            float value = get_value(data[i], j);
                            float yy = up_padding + (1.0f - up_padding - down_padding) * (value - minv[j]) / (maxv[j] - minv[j]);
                            if(!(clicked[j].first <= yy && yy <= clicked[j].second)) is_vaild[i] = 0;
                        }
                    }
                    return true;
                }
            }
            return false;
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
        glm::vec4 change_color(int x, int y, int z, float a)
        {
            return glm::vec4(x / 256.0f, y / 256.0f, z / 256.0f, a);
        }
        void Paint(Common::ImageRGB &input)
        {
            const int width = input.GetSizeX();
            const int height = input.GetSizeY();
            SetBackGround(input, glm::vec4(1.0f, 1.0f, 1.0f, 0.0f));
            ///calculate the maxs and mins
            for (auto &car : data) 
            {
                for(int i = 0; i < (int)(names.size()); ++i)
                {
                    minv[i] = std::min(minv[i], get_value(car, i));
                    maxv[i] = std::max(maxv[i], get_value(car, i));
                }
            }
            std::vector<glm::vec4> colors = {
                change_color(248,128,80,0.4),
                change_color(236,220,10,0.4),
                change_color(129,178,87,0.4),
                change_color(70,225,213,0.4),
                change_color(148,70,225,0.4)
            };
            DrawFilledRect(input, 
                       change_color(133, 245, 10, 0.82),
                        glm::vec2(0.7f, 0.03f),
                        glm::vec2(0.2f, 0.06f)
                    );
            PrintText(input, glm::vec4(0, 0, 0, 1), glm::vec2(0.8, 0.06), 0.03f, "Reset");
            for (size_t i = 0; i < names.size(); ++i) 
            {
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
                if(clicked[i].first != 0 && clicked[i].second != 1)
                {
                    DrawFilledRect(input, 
                        change_color(178,30,208,0.78),
                        glm::vec2(from_x - rect_width / 2, clicked[i].first),
                        glm::vec2(rect_width, clicked[i].second - clicked[i].first)
                    );
                }
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
            int cnt = 0;
            for(int j = 0; j < (int)(data.size()); ++j)
            {
                auto car = data[j];
                glm::vec2 lastPoint;
                cnt++;
                if(!is_vaild[j]) continue;
                for (int i = 0; i < (int)names.size(); ++i) 
                {
                    float x = left_padding + i * line_distance;
                    float value = get_value(car, i);
                    float y = up_padding + (1.0f - up_padding - down_padding) * (value - minv[i]) / (maxv[i] - minv[i]);
                    glm::vec2 point(x, y);
                    if (i > 0) DrawLine(input, colors[cnt % 5], lastPoint, point, 1.0f);
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