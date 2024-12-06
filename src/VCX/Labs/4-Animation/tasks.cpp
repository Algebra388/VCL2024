#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <spdlog/spdlog.h>
#include <iostream>
#include "Labs/4-Animation/tasks.h"
#include "IKSystem.h"
#include "CustomFunc.inl"


namespace VCX::Labs::Animation {
    void ForwardKinematics(IKSystem & ik, int StartIndex) {
        if (StartIndex == 0) {
            ik.JointGlobalRotation[0] = ik.JointLocalRotation[0];
            ik.JointGlobalPosition[0] = ik.JointLocalOffset[0];
            StartIndex                = 1;
        }
        
        for (int i = StartIndex; i < ik.JointLocalOffset.size(); i++) {
            // your code here: forward kinematics, update JointGlobalPosition and JointGlobalRotation
        }
    }

    void InverseKinematicsCCD(IKSystem & ik, const glm::vec3 & EndPosition, int maxCCDIKIteration, float eps) {
        ForwardKinematics(ik, 0);
        // These functions will be useful: glm::normalize, glm::rotation, glm::quat * glm::quat
        for (int CCDIKIteration = 0; CCDIKIteration < maxCCDIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; CCDIKIteration++) {
            // your code here: ccd ik
        }
    }

    void InverseKinematicsFABR(IKSystem & ik, const glm::vec3 & EndPosition, int maxFABRIKIteration, float eps) {
        ForwardKinematics(ik, 0);
        int nJoints = ik.NumJoints();
        std::vector<glm::vec3> backward_positions(nJoints, glm::vec3(0, 0, 0)), forward_positions(nJoints, glm::vec3(0, 0, 0));
        for (int IKIteration = 0; IKIteration < maxFABRIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; IKIteration++) {
            // task: fabr ik
            // backward update
            glm::vec3 next_position         = EndPosition;
            backward_positions[nJoints - 1] = EndPosition;

            for (int i = nJoints - 2; i >= 0; i--) {
                // your code here
            }

            // forward update
            glm::vec3 now_position = ik.JointGlobalPosition[0];
            forward_positions[0] = ik.JointGlobalPosition[0];
            for (int i = 0; i < nJoints - 1; i++) {
                // your code here
            }
            ik.JointGlobalPosition = forward_positions; // copy forward positions to joint_positions
        }

        // Compute joint rotation by position here.
        for (int i = 0; i < nJoints - 1; i++) {
            ik.JointGlobalRotation[i] = glm::rotation(glm::normalize(ik.JointLocalOffset[i + 1]), glm::normalize(ik.JointGlobalPosition[i + 1] - ik.JointGlobalPosition[i]));
        }
        ik.JointLocalRotation[0] = ik.JointGlobalRotation[0];
        for (int i = 1; i < nJoints - 1; i++) {
            ik.JointLocalRotation[i] = glm::inverse(ik.JointGlobalRotation[i - 1]) * ik.JointGlobalRotation[i];
        }
        ForwardKinematics(ik, 0);
    }

    IKSystem::Vec3ArrPtr IKSystem::BuildCustomTargetPosition() {
        // get function from https://www.wolframalpha.com/input/?i=Albert+Einstein+curve
        int nums = 5000;
        using Vec3Arr = std::vector<glm::vec3>;
        std::shared_ptr<Vec3Arr> custom(new Vec3Arr(nums));
        int index = 0;
        for (int i = 0; i < nums; i++) {
            float x_val = 1.5e-3f * custom_x(92 * glm::pi<float>() * i / nums);
            float y_val = 1.5e-3f * custom_y(92 * glm::pi<float>() * i / nums);
            if (std::abs(x_val) < 1e-3 || std::abs(y_val) < 1e-3) continue;
            (*custom)[index++] = glm::vec3(1.6f - x_val, 0.0f, y_val - 0.2f);
        }
        custom->resize(index);
        return custom;
    }

    static Eigen::VectorXf glm2eigen(std::vector<glm::vec3> const & glm_v) {
        Eigen::VectorXf v = Eigen::Map<Eigen::VectorXf const, Eigen::Aligned>(reinterpret_cast<float const *>(glm_v.data()), static_cast<int>(glm_v.size() * 3));
        return v;
    }

    static std::vector<glm::vec3> eigen2glm(Eigen::VectorXf const & eigen_v) {
        return std::vector<glm::vec3>(
            reinterpret_cast<glm::vec3 const *>(eigen_v.data()),
            reinterpret_cast<glm::vec3 const *>(eigen_v.data() + eigen_v.size())
        );
    }

    static Eigen::SparseMatrix<float> CreateEigenSparseMatrix(std::size_t n, std::vector<Eigen::Triplet<float>> const & triplets) {
        Eigen::SparseMatrix<float> matLinearized(n, n);
        matLinearized.setFromTriplets(triplets.begin(), triplets.end());
        return matLinearized;
    }

    // solve Ax = b and return x
    static Eigen::VectorXf ComputeSimplicialLLT(
        Eigen::SparseMatrix<float> const & A,
        Eigen::VectorXf const & b) {
        auto solver = Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>(A);
        return solver.solve(b);
    }

    void AdvanceMassSpringSystem(MassSpringSystem & system, float const dt) {
        // your code here: rewrite following code
        int const steps = 1000;
        float const ddt = dt / steps; 
        for (std::size_t s = 0; s < steps; s++) {
            /* explicit
            std::vector<glm::vec3> forces(system.Positions.size(), glm::vec3(0));
            for (auto const spring : system.Springs) { /// every spring
                auto const p0 = spring.AdjIdx.first;
                auto const p1 = spring.AdjIdx.second;
                glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
                glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
                glm::vec3 const e01 = glm::normalize(x01);
                glm::vec3 f = (system.Stiffness * (glm::length(x01) - spring.RestLength) + system.Damping * glm::dot(v01, e01)) * e01;
                forces[p0] += f;
                forces[p1] -= f;
            }
            for (std::size_t i = 0; i < system.Positions.size(); i++) {
                if (system.Fixed[i]) continue;
                system.Velocities[i] += (glm::vec3(0, -system.Gravity, 0) + forces[i] / system.Mass) * ddt;
                system.Positions[i] += system.Velocities[i] * ddt;
            }

            */

           /* implicit */
           std::vector<glm::vec3> matrix_g(system.Positions.size(), glm::vec3(0));
           std::vector<glm::mat3> matrix_hg_accumulation_dignal(system.Positions.size(), glm::mat3(0));
           std::vector<Eigen::Triplet<float>> triples;/// for HGs
           std::vector<glm::vec3> f_ext(system.Positions.size(), glm::vec3(0));
           /// this is optional for parallel
           const glm::mat3 M {system.Mass, .0, .0, .0, system.Mass, .0, .0, .0, system.Mass};
           const glm::mat3 Indentity {1.0, .0, .0, .0, 1.0, .0, .0, .0, 1.0};
           /// calculation of g
           /// calculation of hg
           for (auto const spring : system.Springs) { /// every spring
                auto const p0 = spring.AdjIdx.first;
                auto const p1 = spring.AdjIdx.second;
                auto const pos0 = system.Positions[p0];
                auto const pos1 = system.Positions[p1];
                glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
                glm::vec3 const e01 = glm::normalize(x01);
                glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
                const float length = glm::length(x01);
                const float damping = system.Damping; /// external force
                const float origional_length = spring.RestLength;
                const float kij = system.Stiffness;
                /// calculate external force
                f_ext[p0] += system.Damping * glm::dot(v01, e01) * e01;
                f_ext[p1] -= system.Damping * glm::dot(v01, e01) * e01;
                ///calculate E for g
                ///carefule +-
                matrix_g[p0] -= kij * (length - origional_length) * e01;
                matrix_g[p1] += kij * (length - origional_length) * e01;
                glm::mat3 H_current {.0f};
                H_current = kij * glm::outerProduct(x01, x01) / (length * length) +\
                    kij * (1 - origional_length / length) * (Indentity - \
                    glm::outerProduct(x01, x01) / (length * length));
                matrix_hg_accumulation_dignal[p0] += H_current;
                matrix_hg_accumulation_dignal[p1] += H_current;
                for(int i = 0; i < 3; ++i)
                {
                    for(int j = 0; j < 3; ++j)
                    {/// the elements not in diag
                        triples.push_back({3 * p0 + i, 3 * p1 + j, - H_current[i][j]});
                        triples.push_back({3 * p1 + i, 3 * p0 + j, - H_current[i][j]});
                    }
                }
            }
           for(std::size_t i = 0; i < system.Positions.size(); ++i)
           {
                // part1 : HG
                if (system.Fixed[i]) 
                {
                    for (int j = 0; j < 3; ++j) 
                    {
                        triples.push_back({i * 3 + j, i * 3 + j, system.Mass / (ddt * ddt) });
                    }
                    continue;
                }/// not move only additional
                for(int p = 0; p < 3; ++p)
                {
                    for(int q = 0; q < 3; ++q)
                    {
                        float additional_value = 0;
                        if(p == q) additional_value = system.Mass / ddt / ddt;
                        triples.push_back({i * 3 + p, i * 3 + q, matrix_hg_accumulation_dignal[i][p][q] + additional_value});
                    }
                }
                glm::vec3 yki = system.Positions[i] + system.Velocities[i] * ddt + ddt * \
                ddt / system.Mass * f_ext[i] + glm::vec3(0, -system.Gravity, 0);
                matrix_g[i] += system.Mass / ddt / ddt * (system.Positions[i] - yki);
           }
           auto A = CreateEigenSparseMatrix(system.Positions.size() * 3, triples);
           auto b = glm2eigen(matrix_g);
           auto x = ComputeSimplicialLLT(A, b);
           std::vector<glm::vec3> delta_x = eigen2glm(x);
           for(std::size_t i = 0; i < system.Positions.size(); ++i)
           {
                if (system.Fixed[i]) 
                {
                    system.Velocities[i] = glm::vec3(0);
                    continue;
                }
                system.Velocities[i] = delta_x[i] / ddt;
                system.Positions[i] += delta_x[i];
           }
        }
    }
}
