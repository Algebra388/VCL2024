#include "Engine/app.h"
#include "Labs/4-Animation/CaseMassSpring.h"
#include "Labs/4-Animation/tasks.h"
#include "Labs/Common/ImGuiHelper.h"
#include <random>
#include <map>

namespace VCX::Labs::Animation {
    CaseMassSpring::CaseMassSpring() :
        _program(
            Engine::GL::UniqueProgram({
                Engine::GL::SharedShader("assets/shaders/flat.vert"),
                Engine::GL::SharedShader("assets/shaders/flat.frag") })),
        _particlesItem(Engine::GL::VertexLayout()
            .Add<glm::vec3>("position", Engine::GL::DrawFrequency::Stream , 0), Engine::GL::PrimitiveType::Points),
        _springsItem(Engine::GL::VertexLayout()
            .Add<glm::vec3>("position", Engine::GL::DrawFrequency::Stream , 0), Engine::GL::PrimitiveType::Lines) {
        _cameraManager.AutoRotate = false;
        _cameraManager.Save(_camera);
        ResetSystem();
    }

    void CaseMassSpring::OnSetupPropsUI() {
        if (ImGui::CollapsingHeader("Model", ImGuiTreeNodeFlags_DefaultOpen)) {
            const char * model_option[] = { "square fabric", "trampoline"};
            ///bool pressed = 0;
            if (ImGui::BeginCombo("Model op", model_option[_massSpringSystem.model_type])) 
            {
                for (int i = 0; i < 2; i++) {
                    bool selected = i == _massSpringSystem.model_type;
                    if (ImGui::Selectable(model_option[i], selected)) {
                        if (! selected) 
                        {
                            _massSpringSystem.model_type = i;
                            Reload();
                        }
                    }
                }
                ImGui::EndCombo();
            }
        }
        if (ImGui::CollapsingHeader("Algorithm", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ImGui::Button("Reset System")) ResetSystem();
            ImGui::SameLine();
            if (ImGui::Button(_stopped ? "Start Simulation" : "Stop Simulation")) _stopped = ! _stopped;
            const char * Mass_options[] = { "uniform distribution", "same mass"};
            ///bool pressed = 0;
            if (ImGui::BeginCombo("Mass options", Mass_options[_massSpringSystem.Mass_type])) 
            {
                for (int i = 0; i < 2; i++) {
                    bool selected = i == _massSpringSystem.Mass_type;
                    if (ImGui::Selectable(Mass_options[i], selected)) {
                        if (! selected) 
                        {
                            _massSpringSystem.Mass_type = i;
                        }
                    }
                }
                ImGui::EndCombo();
            }
            if( _massSpringSystem.Mass_type == 1)
            {
                ImGui::SliderFloat("Part. Mass", &_massSpringSystem.mass, .5f, 10.f);
                _massSpringSystem.lower = _massSpringSystem.upper = _massSpringSystem.mass;
            }
            else
            {
                ImGui::SliderFloat("Part. Mass Lower Bound", &_massSpringSystem.lower, .5f, 10.f);
                ImGui::SliderFloat("Part. Mass Upper Bound", &_massSpringSystem.upper, .5f, 10.f);
            }
            std::uniform_real_distribution<> dist(_massSpringSystem.lower, _massSpringSystem.upper);
            std::mt19937 rng(20250115);
            /// use same seed, if the parameter is same, got the same mass each node
            _massSpringSystem.Mass.clear();
            for (std::size_t i = 0; i < _massSpringSystem.vertexs; i++) {
                if(_massSpringSystem.Mass_type == 0)
                {
                    _massSpringSystem.Mass.push_back(dist(rng));
                }
                else _massSpringSystem.Mass.push_back(1);
            }
            //printf("%f\n",_massSpringSystem.Mass[0]);
            ImGui::SliderInt("Times of Iterations", &_massSpringSystem.Iterations, 1, 100);
            ImGui::SliderFloat("External Force X-axis", &_massSpringSystem.Ext_force.x, -10.f, 10.f);
            ImGui::SliderFloat("External Force Y-axis", &_massSpringSystem.Ext_force.y, -10.f, 10.f);
            ImGui::SliderFloat("External Force Z-axis", &_massSpringSystem.Ext_force.z, -10.f, 10.f);
            //ImGui::SliderFloat("Spr. Damp.", &_massSpringSystem.Damping, .1f, 10.f);
            ImGui::SliderFloat("Gravity", &_massSpringSystem.Gravity, .1f, 1.f);
        }
        ImGui::Spacing();
        if (ImGui::CollapsingHeader("Appearance")) {
            ImGui::SliderFloat("Part. Size", &_particleSize, 1, 6);
            ImGui::ColorEdit3("Part. Color", glm::value_ptr(_particleColor));
            ImGui::SliderFloat("Spr. Width", &_springWidth, .001f, 1.f);
            ImGui::ColorEdit3("Spr. Color", glm::value_ptr(_springColor));
        }
        ImGui::Spacing();
    }

    Common::CaseRenderResult CaseMassSpring::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        if (! _stopped) AdvanceMassSpringSystem(_massSpringSystem, Engine::GetDeltaTime());
        
        _particlesItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(_massSpringSystem.Positions));
        _springsItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(_massSpringSystem.Positions));

        _frame.Resize(desiredSize);

        _cameraManager.Update(_camera);

        _program.GetUniforms().SetByName("u_Projection", _camera.GetProjectionMatrix((float(desiredSize.first) / desiredSize.second)));
        _program.GetUniforms().SetByName("u_View"      , _camera.GetViewMatrix());

        gl_using(_frame);
        glEnable(GL_LINE_SMOOTH);
        glPointSize(_particleSize);
        glLineWidth(_springWidth);

        _program.GetUniforms().SetByName("u_Color", _springColor);
        _springsItem.Draw({ _program.Use() });
        _program.GetUniforms().SetByName("u_Color", _particleColor);
        _particlesItem.Draw({ _program.Use() });

        glLineWidth(1.f);
        glPointSize(1.f);
        glDisable(GL_LINE_SMOOTH);

        return Common::CaseRenderResult {
            .Fixed     = false,
            .Flipped   = true,
            .Image     = _frame.GetColorAttachment(),
            .ImageSize = desiredSize,
        };
    }
    void CaseMassSpring::OnProcessInput(ImVec2 const & pos) {
        _cameraManager.ProcessInput(_camera, pos);
    }

    void CaseMassSpring::Reload()
    {
        int type = _massSpringSystem.model_type;
        _massSpringSystem = { };
        _massSpringSystem.model_type = type;
        //puts("1");
        if(_massSpringSystem.model_type == 0)
        {
            std::size_t const n = 10;
            float const delta = 2.f / n;
            _massSpringSystem.vertexs = (n + 1) * (n + 1); 
            auto constexpr GetID = [](std::size_t const i, std::size_t const j) { return i * (n + 1) + j; };
            for (std::size_t i = 0; i <= n; i++) {
                for (std::size_t j = 0; j <= n; j++) {
                    _massSpringSystem.AddParticle(glm::vec3(i * delta , 1.5f, j * delta - 1.f));
                    if (i > 0) _massSpringSystem.AddSpring(GetID(i, j), GetID(i - 1, j));
                    if (i > 1) _massSpringSystem.AddSpring(GetID(i, j), GetID(i - 2, j));
                    if (j > 0) _massSpringSystem.AddSpring(GetID(i, j), GetID(i, j - 1));
                    if (j > 1) _massSpringSystem.AddSpring(GetID(i, j), GetID(i, j - 2));
                    if (i > 0 && j > 0) _massSpringSystem.AddSpring(GetID(i, j), GetID(i - 1, j - 1));
                    if (i > 0 && j < n) _massSpringSystem.AddSpring(GetID(i, j), GetID(i - 1, j + 1));
                }
            }
            _massSpringSystem.Fixed[GetID(0, 0)] = true;
            _massSpringSystem.Fixed[GetID(0, n)] = true;
            // _massSpringSystem.Fixed[GetID(n, 0)] = true;
            // _massSpringSystem.Fixed[GetID(n, n)] = true;
            
        }
        else /// a circle(maybe)
        {
            //puts("1");
            std::size_t const out_side = 20;/// divisible by 4
            int radis = 4;
            //std::map<std::pair<float, float>, int>index;
            _massSpringSystem.mass = 1;
            _massSpringSystem.upper = _massSpringSystem.lower = 1;
            std::vector<glm::vec3> position;
            float z = 1.5;
            for(int i = 0; i < out_side; ++i)
            {
                float x = glm::cos(2 * glm::acos(-1) / out_side * i) * radis;
                float y = glm::sin(2 * glm::acos(-1) / out_side * i) * radis;
                //index[{x, y}] = i;
                _massSpringSystem.AddParticle(glm::vec3(x, z, y));
                position.push_back(glm::vec3(x, z, y));
                _massSpringSystem.Fixed[i] = true;
                //printf("%d = (%f, %f)\n", i, x, y);
            }
            for(int i = 0; i < out_side; ++i) _massSpringSystem.AddSpring(i, (i + 1) % out_side);
            //puts("2");
            int tot = out_side;
            for(int i = 1; i < out_side / 2; ++i)
            {
                for(int j = out_side / 4 + 1; j < out_side / 4 * 3; ++j)
                {
                    //index[{position[j].x, positions[i].y}] = tot++;
                    auto t = glm::vec3(position[i].x, z, position[j].z);
                    if(glm::length(t) > radis) continue;
                    position.push_back(t);
                    t.y = z - glm::sqrt(radis * radis - glm::length(t) * glm::length(t));
                    _massSpringSystem.AddParticle(t);
                    //printf("%f %f\n",position[i].x, position[j].z);
                    tot++;
                }
            }
            //puts("3");
            for(int i = 0; i < tot; ++i)
            {
                float x = position[i].x, y = position[i].z;
                int idx = -1, idy = -1;
                for(int j = 0; j < tot; ++j)
                {
                    float nx = position[j].x, ny = position[j].z;
                    if(nx == x && ny > y)
                    {
                        if(idx == -1 || glm::length(position[i] - position[j]) < glm::length(position[i] - position[idx]))
                        {
                            idx = j;
                        }
                    }
                    if(nx > x && ny == y)
                    {
                        if(idy == -1 || glm::length(position[i] - position[j]) < glm::length(position[i] - position[idy]))
                        {
                            idy = j;
                        }
                    }
                }
                if(idx != -1) _massSpringSystem.AddSpring(i, idx); //printf("%d %d\n", i, idx);
                if(idy != -1) _massSpringSystem.AddSpring(i, idy); //printf("%d %d\n", i, idy);
            }
            //puts("4");
            _massSpringSystem.vertexs = tot;
        }
        
        _massSpringSystem.mass = 1;
        _massSpringSystem.upper = _massSpringSystem.lower = 1;
        for (std::size_t i = 0; i <= _massSpringSystem.vertexs; i++) {
            _massSpringSystem.Mass.push_back(_massSpringSystem.mass);
        }
        std::vector<std::uint32_t> indices;
        for (auto const & spring : _massSpringSystem.Springs) {
            indices.push_back(std::uint32_t(spring.AdjIdx.first));
            indices.push_back(std::uint32_t(spring.AdjIdx.second));
        }
        _springsItem.UpdateElementBuffer(indices);
    }
    void CaseMassSpring::ResetSystem() {
        _massSpringSystem.model_type = 0;
        Reload();
    }
}
