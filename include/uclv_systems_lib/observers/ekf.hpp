/*
    System interface Class Discrete Time System

    Copyright 2024 Università della Campania Luigi Vanvitelli

    Authors: Marco Costanzo  <marco.costanzo@unicampania.it>
             Marco De Simone <marco.desimone@unicampania.it>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <memory>
#include <Eigen/Dense>
#include "../ss/state_space_interface.hpp"

namespace uclv::systems
{

  template <typename Scalar_t, int dim_state, int dim_input, int dim_output>
  class ExtendedKalmanFilter
  {
  public:
    using SharedPtr = std::shared_ptr<ExtendedKalmanFilter>;
    using ConstSharedPtr = std::shared_ptr<const ExtendedKalmanFilter>;
    using WeakPtr = std::weak_ptr<ExtendedKalmanFilter>;
    using ConstWeakPtr = std::weak_ptr<const ExtendedKalmanFilter>;
    using UniquePtr = std::unique_ptr<ExtendedKalmanFilter>;

    using StateSpaceInterface_t = StateSpaceInterface<Scalar_t, dim_state, dim_input, dim_output>;

    using Input_t = typename StateSpaceInterface_t::Input_t;
    using InputRef_t = typename StateSpaceInterface_t::InputRef_t;
    using InputConstRef_t = typename StateSpaceInterface_t::InputConstRef_t;
    using Output_t = typename StateSpaceInterface_t::Output_t;
    using OutputRef_t = typename StateSpaceInterface_t::OutputRef_t;
    using OutputConstRef_t = typename StateSpaceInterface_t::OutputConstRef_t;

    using State_t = typename StateSpaceInterface_t::State_t;
    using StateRef_t = typename StateSpaceInterface_t::StateRef_t;
    using StateConstRef_t = typename StateSpaceInterface_t::StateConstRef_t;

    using JacobianStateState_t = typename StateSpaceInterface_t::JacobianStateState_t;
    using JacobianStateStateRef_t = typename StateSpaceInterface_t::JacobianStateStateRef_t;
    using JacobianStateStateConstRef_t = typename StateSpaceInterface_t::JacobianStateStateConstRef_t;

    using JacobianStateInput_t = typename StateSpaceInterface_t::JacobianStateInput_t;
    using JacobianStateInputRef_t = typename StateSpaceInterface_t::JacobianStateInputRef_t;
    using JacobianStateInputConstRef_t = typename StateSpaceInterface_t::JacobianStateInputConstRef_t;

    using JacobianOutputState_t = typename StateSpaceInterface_t::JacobianOutputState_t;
    using JacobianOutputStateRef_t = typename StateSpaceInterface_t::JacobianOutputStateRef_t;
    using JacobianOutputStateConstRef_t = typename StateSpaceInterface_t::JacobianOutputStateConstRef_t;

    using JacobianOutputInput_t = typename StateSpaceInterface_t::JacobianOutputInput_t;
    using JacobianOutputInputRef_t = typename StateSpaceInterface_t::JacobianOutputInputRef_t;
    using JacobianOutputInputConstRef_t = typename StateSpaceInterface_t::JacobianOutputInputConstRef_t;

    using StateNoiseCovariance_t = Eigen::Matrix<Scalar_t, dim_state, dim_state>;
    using StateNoiseCovarianceRef_t = Eigen::Ref<StateNoiseCovariance_t>;
    using StateNoiseCovarianceConstRef_t = Eigen::Ref<const StateNoiseCovariance_t>;

    using StateCovariance_t = StateNoiseCovariance_t;
    using StateCovarianceRef_t = Eigen::Ref<StateCovariance_t>;
    using StateCovarianceConstRef_t = Eigen::Ref<const StateCovariance_t>;

    using OutputNoiseCovariance_t = Eigen::Matrix<Scalar_t, dim_output, dim_output>;
    using OutputNoiseCovarianceRef_t = Eigen::Ref<OutputNoiseCovariance_t>;
    using OutputNoiseCovarianceConstRef_t = Eigen::Ref<const OutputNoiseCovariance_t>;

    using NormalizationFunction_t = void (*)(const StateConstRef_t &, StateRef_t);

    ExtendedKalmanFilter(typename StateSpaceInterface_t::SharedPtr system_ptr,
                         const StateNoiseCovarianceConstRef_t &W,
                         const OutputNoiseCovarianceConstRef_t &V,
                         NormalizationFunction_t normalization_fun = nullptr)
        : system_(system_ptr), W_(W), V_(V), P_(W),
          normalization_fun_(normalization_fun)
    {
      Identity_x_.resizeLike(P_);
      Identity_x_.setIdentity();
      x_hat_k_k_.resizeLike(system_ptr->get_state());
      y_hat_k_.resizeLike(system_ptr->get_output());
    }

    ExtendedKalmanFilter(const ExtendedKalmanFilter &other)
        : system_(other.system_->clone()), P_(other.P_), W_(other.W_), V_(other.V_), Identity_x_(other.Identity_x_), x_hat_k_k_(other.x_hat_k_k_), y_hat_k_(other.y_hat_k_)
    {
    }

    void setP(const StateCovarianceConstRef_t &P)
    {
      P_ = P;
    }

    void setW(const StateNoiseCovarianceConstRef_t &W)
    {
      W_ = W;
    }

    void setV(const OutputNoiseCovarianceConstRef_t &V)
    {
      V_ = V;
    }

    void kf_apply(const InputConstRef_t &u_k,
                  const OutputConstRef_t &y_k,
                  const StateNoiseCovarianceConstRef_t &W_k,
                  const OutputNoiseCovarianceConstRef_t &V_k)
    {
      setW(W_k);
      setV(V_k);
      obs_apply(u_k, y_k);
    }

    void obs_apply(const InputConstRef_t &u_k,
                   const OutputConstRef_t &y_k)
    {
      // PREDICT
      State_t x_hat_k_k1;
      system_->state_fcn(x_hat_k_k_, u_k, x_hat_k_k1);

      if (normalization_fun_ != nullptr)
      {
        normalization_fun_(x_hat_k_k1, x_hat_k_k1);
      }

      JacobianStateState_t F_k1;
      system_->jacobx_state_fcn(x_hat_k_k_, u_k, F_k1);
      StateCovariance_t P_k_k1 = F_k1 * P_ * F_k1.transpose() + W_;

      // UPDATE
      Output_t y_hat_k_k1;
      system_->output_fcn(x_hat_k_k1, u_k, y_hat_k_k1);
      Output_t y_tilde_k = y_k - y_hat_k_k1;

      JacobianOutputState_t H_k;
      system_->jacobx_output_fcn(x_hat_k_k1, u_k, H_k);
      OutputNoiseCovariance_t S_k = H_k * P_k_k1 * H_k.transpose() + V_;
      auto K_k = P_k_k1 * H_k.transpose() * S_k.inverse();
      x_hat_k_k_ = x_hat_k_k1 + K_k * y_tilde_k;

      if (normalization_fun_ != nullptr)
      {
        normalization_fun_(x_hat_k_k_, x_hat_k_k_);
      }
      P_ = (Identity_x_ - K_k * H_k) * P_k_k1;
      system_->output_fcn(x_hat_k_k_, u_k, y_hat_k_);
    }

    const State_t &get_state() const
    {
      return x_hat_k_k_;
    }
    const Output_t &get_output() const
    {
      return y_hat_k_;
    }

    void set_state(const StateRef_t &x)
    {
      x_hat_k_k_ = x;
    }

    void reset()
    {
      P_ = W_;
      system_->get_resetted_state(x_hat_k_k_);
      system_->get_resetted_output(y_hat_k_);
    }

    void display()
    {
      std::cout << "\n"
                << std::endl;
      std::cout << "Extended Kalman Filter\n"
                << std::endl;
      std::cout << "W: \n"
                << W_ << std::endl;
      std::cout << "V: \n"
                << V_ << std::endl;
      system_->display();
    }

  private:
    typename StateSpaceInterface_t::SharedPtr system_;
    StateNoiseCovariance_t W_;
    OutputNoiseCovariance_t V_;
    StateCovariance_t P_;
    StateCovariance_t Identity_x_;
    State_t x_hat_k_k_;
    Output_t y_hat_k_;
    NormalizationFunction_t normalization_fun_ = nullptr;
  };

} // namespace uclv::systems