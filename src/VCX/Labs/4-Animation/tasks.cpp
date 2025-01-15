#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <spdlog/spdlog.h>
#include <iostream>
#include "Labs/4-Animation/tasks.h"
#include "IKSystem.h"
#include "CustomFunc.inl"


namespace VCX::Labs::Animation {

   /*
  std::vector<glm::vec3>   JointLocalOffset;
        std::vector<float>       JointOffsetLength;
        std::vector<glm::vec3>   JointGlobalPosition;
        std::vector<glm::quat>   JointLocalRotation;
        std::vector<glm::quat>   JointGlobalRotation; 
   */

    /*
     实现前向运动学 (forward kinematic) 算法。程序入口在 tasks.cpp 中 
     ForwardKinematics 函数，参数 StartIndex 给出了前向更新从第几个关节开始的信息。 
     04 - A: 函数 glm::rotate(q, v) 或者直接使用重载后的乘号 q * v 即为向量 v 经过四元数 q 旋转后的结果。 
    */
    void ForwardKinematics(IKSystem & ik, int StartIndex) {
        if (StartIndex == 0) {
            ik.JointGlobalRotation[0] = ik.JointLocalRotation[0];
            ik.JointGlobalPosition[0] = ik.JointLocalOffset[0];
            StartIndex                = 1;
        }
         /// using quad number
        for (int i = StartIndex; i < ik.JointLocalOffset.size(); i++) {
            // your code here: forward kinematics, update JointGlobalPosition and JointGlobalRotation
            ik.JointGlobalRotation[i] = ik.JointGlobalRotation[i - 1] * ik.JointLocalRotation[i];
            ik.JointGlobalPosition[i] = ik.JointGlobalPosition[i - 1] + \
                ik.JointGlobalRotation[i - 1] * ik.JointLocalOffset[i];
        }
    }

    void InverseKinematicsCCD(IKSystem & ik, const glm::vec3 & EndPosition, int maxCCDIKIteration, float eps) {
        ForwardKinematics(ik, 0);
        int CCDIKIteration;
        // These functions will be useful: glm::normalize, glm::rotation, glm::quat * glm::quat
        for (CCDIKIteration = 0; CCDIKIteration < maxCCDIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; CCDIKIteration++) {
            // your code here: ccd ik
            int n = ik.JointLocalOffset.size() - 1;/// points from 0 to n
            for(int i = n - 1; i >= 0; --i)
            {
                glm::vec3 endp = ik.JointGlobalPosition[n];
                glm::vec3 startp = ik.JointGlobalPosition[i];
                glm::vec3 target = glm::normalize(EndPosition - startp);
                glm::vec3 current = glm::normalize(endp - startp);
                auto angle = glm::rotation(current, target);
                ik.JointLocalRotation[i] = angle * ik.JointLocalRotation[i];
                ForwardKinematics(ik, i);
            }
        }
        //printf("%d\n", CCDIKIteration);
    }

    void InverseKinematicsFABR(IKSystem & ik, const glm::vec3 & EndPosition, int maxFABRIKIteration, float eps) {
        ForwardKinematics(ik, 0);
        int nJoints = ik.NumJoints();
        int IKIteration;
        std::vector<glm::vec3> backward_positions(nJoints, glm::vec3(0, 0, 0)), forward_positions(nJoints, glm::vec3(0, 0, 0));
        for (IKIteration = 0; IKIteration < maxFABRIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; IKIteration++) {
            // task: fabr ik
            // backward update
            glm::vec3 next_position         = EndPosition;
            backward_positions[nJoints - 1] = EndPosition;

            for (int i = nJoints - 2; i >= 0; i--) {
                // your code here
                glm::vec3 dirc = glm::normalize(next_position - ik.JointGlobalPosition[i]);
                glm::vec3 nowp = next_position - dirc * ik.JointOffsetLength[i + 1];
                backward_positions[i] = nowp;
                next_position = backward_positions[i];
            }

            // forward update
            glm::vec3 now_position = ik.JointGlobalPosition[0];
            forward_positions[0] = ik.JointGlobalPosition[0];
            for (int i = 0; i < nJoints - 1; i++) {
                // your code here
                glm::vec3 dirc = glm::normalize(backward_positions[i + 1] - now_position);
                forward_positions[i + 1] = now_position + dirc * ik.JointOffsetLength[i + 1];
                now_position = forward_positions[i + 1];
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
        //printf("%d\n", IKIteration);
    }

    float custom_xx(float t)
    {
        if(t >= 0 && t <= 0.1) return t * 2;
        else if(t >= 0.4 && t <= 0.6) return (0.5 - t) * 2;
        else if(t >= 0.1 && t <= 0.4) return 0.2;
        else if(t >= 0.6 && t <= 0.9) return -0.2;
        else return (t - 1) * 2;
    }

    float custom_yy(float t)
    {
        if(t >= 0 && t <= 0.1) return -0.3;
        else if(t >= 0.4 && t <= 0.6) return 0.3;
        else if(t >= 0.1 && t <= 0.4) return (t - 0.25) * 2;
        else if(t >= 0.6 && t <= 0.9) return (0.75 - t) * 2;
        else return -0.3;
    }

    IKSystem::Vec3ArrPtr IKSystem::BuildCustomTargetPosition() {
        // get function from https://www.wolframalpha.com/input/?i=Albert+Einstein+curve
        int nums = 100;
        using Vec3Arr = std::vector<glm::vec3>;
        std::shared_ptr<Vec3Arr> custom(new Vec3Arr(nums));
        int index = 0;
        /*for (int i = 0; i < nums; i++) {
            float x_val = 1.5e-3f * custom_xx(92 * glm::pi<float>() * i / nums);
            float y_val = 1.5e-3f * custom_yy(92 * glm::pi<float>() * i / nums);
            if (std::abs(x_val) < 1e-3 || std::abs(y_val) < 1e-3) continue;
            (*custom)[index++] = glm::vec3(1.6f - x_val, 0.0f, y_val - 0.2f);
        }*/
       for(int i = 0; i < nums; ++i)
       {
            float x_val = custom_xx(1.0 * i / nums);
            float y_val = custom_yy(1.0 * i / nums);
            if (std::abs(x_val) < 1e-3 || std::abs(y_val) < 1e-3) continue;
            (*custom)[index++] = glm::vec3(x_val, 0.0f, y_val);
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
        //int const steps = 1000;
        int const steps = 10; /// times of PBD
        float const ddt = dt / steps; 
        double st = clock();
        for (std::size_t s = 0; s < steps; s++) {
            //explicit
            std::vector<glm::vec3> temp_Positions(system.Positions.size(), glm::vec3(0));
            /*for (auto const spring : system.Springs) { /// every spring
                auto const p0 = spring.AdjIdx.first;
                auto const p1 = spring.AdjIdx.second;
                glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
                glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
                glm::vec3 const e01 = glm::normalize(x01);
                glm::vec3 f = (system.Stiffness * (glm::length(x01) - spring.RestLength) + system.Damping * glm::dot(v01, e01)) * e01;
                forces[p0] += f;
                forces[p1] -= f;
            }*/
            for (std::size_t i = 0; i < system.Positions.size(); i++) {
                if (system.Fixed[i]) 
                {
                    temp_Positions[i] = system.Positions[i];
                    continue;
                }
                //system.Velocities[i] += (glm::vec3(0, -system.Gravity, 0) + forces[i] / system.Mass) * ddt;
                temp_Positions[i] = system.Positions[i] + system.Velocities[i] * ddt;
                temp_Positions[i] += (glm::vec3(0, -system.Gravity, 0) + system.Ext_force)* ddt * ddt / 2.0f;
            }//external gravity

            ///Iterations

            for(int iter = 0; iter < system.Iterations; ++iter)
            {
                for(auto const &spring : system.Springs)///for every spring
                {
                    auto const index1 = spring.AdjIdx.first;
                    auto const index2 = spring.AdjIdx.second;/// the index
                    auto const p1 = temp_Positions[index1];
                    auto const p2 = temp_Positions[index2];
                    auto const mass1 = system.Mass[index1];
                    auto const mass2 = system.Mass[index2];
                    //printf("%f %f\n", mass1, mass2);
                    auto const length = spring.RestLength;
                    auto delta_p1 = - (mass2 / (mass1 + mass2)) * (glm::length(p1 - p2) - length) * (p1 - p2) / glm::length(p1 - p2);
                    auto delta_p2 = mass1 / (mass1 + mass2) * (glm::length(p1 - p2) - length) * (p1 - p2) / glm::length(p1 - p2);
                    if(!system.Fixed[index1])temp_Positions[index1] += delta_p1;
                    if(!system.Fixed[index2])temp_Positions[index2] += delta_p2;
                }
            }
            for (std::size_t i = 0; i < system.Positions.size(); i++)
            {
                if(system.Fixed[i])
                {
                    system.Velocities[i] = glm::vec3(0);
                    continue;
                }
                system.Velocities[i] = (temp_Positions[i] - system.Positions[i]) / ddt;
                system.Positions[i] = temp_Positions[i];
            }
            /*
           printf("Step %d\n", s);
           for (std::size_t i = 0; i < system.Positions.size(); i++) {
                printf("Position[%d] = %.2lf %.2lf %.2lf\n", i, system.Positions[i].x, system.Positions[i].y, system.Positions[i].z);
                printf("Velocities[%d] = %.2lf %.2lf %.2lf\n", i, system.Velocities[i].x, system.Velocities[i].y, system.Velocities[i].z);
                //system.Velocities[i] 
                //system.Positions[i] += system.Velocities[i] * ddt;
            }
            */
        }
        //while(clock() - st <= 5000);
    }
}