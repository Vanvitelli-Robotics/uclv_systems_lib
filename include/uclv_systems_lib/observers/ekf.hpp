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

template <int dim_state, int dim_input, int dim_output>
class ExtendedKalmanFilter
{
public:
  typedef std::shared_ptr<ExtendedKalmanFilter> SharedPtr;
  typedef std::shared_ptr<const ExtendedKalmanFilter> ConstSharedPtr;
  typedef std::weak_ptr<ExtendedKalmanFilter> WeakPtr;
  typedef std::weak_ptr<const ExtendedKalmanFilter> ConstWeakPtr;
  typedef std::unique_ptr<ExtendedKalmanFilter> UniquePtr;

  typedef ::uclv::systems::StateSpaceInterface<dim_state, dim_input, dim_output> StateSpaceInterface;

  ExtendedKalmanFilter(typename StateSpaceInterface::SharedPtr system_ptr,
                       const Eigen::Matrix<double, dim_state, dim_state>& W,
                       const Eigen::Matrix<double, dim_output, dim_output>& V)
    : system_(system_ptr), W_(W), V_(V), P_(W)
  {
    Identity_x_.resizeLike(P_);
    Identity_x_.setIdentity();
    x_hat_k_k_.resizeLike(system_ptr->get_state());
    y_hat_k_.resizeLike(system_ptr->get_output());
  }

  ExtendedKalmanFilter(const ExtendedKalmanFilter& other)
    : system_(other.system_->clone())
    , P_(other.P_)
    , W_(other.W_)
    , V_(other.V_)
    , Identity_x_(other.Identity_x_)
    , x_hat_k_k_(other.x_hat_k_k_)
    , y_hat_k_(other.y_hat_k_)
  {
  }

  void setP(const Eigen::Matrix<double, dim_state, dim_state>& P)
  {
    P_ = P;
  }

  void setW(const Eigen::Matrix<double, dim_state, dim_state>& W)
  {
    W_ = W;
  }

  void setV(const Eigen::Matrix<double, dim_output, dim_output>& V)
  {
    V_ = V;
  }

  void kf_apply(const Eigen::Ref<const Eigen::Matrix<double, dim_input, 1>>& u_k,
                const Eigen::Ref<const Eigen::Matrix<double, dim_output, 1>>& y_k,
                const Eigen::Matrix<double, dim_state, dim_state>& W_k,
                const Eigen::Matrix<double, dim_output, dim_output>& V_k)
  {
    setW(W_k);
    setV(V_k);
    obs_apply(u_k, y_k);
  }

  void obs_apply(const Eigen::Ref<const Eigen::Matrix<double, dim_input, 1>>& u_k,
                 const Eigen::Ref<const Eigen::Matrix<double, dim_output, 1>>& y_k)
  {
    // Eigen::Matrix<double, dim_state, 1> x_hat_k1_k1 = x_hat_k_k_;
    Eigen::Matrix<double, dim_state, dim_state> P_k1_k1 = P_;
    Eigen::Matrix<double, dim_state, dim_state> W_k1 = W_;
    Eigen::Matrix<double, dim_output, dim_output> V_k = V_;

    // PREDICT
    Eigen::Matrix<double, dim_state, 1> x_hat_k_k1;
    system_->state_fcn(x_hat_k_k_, u_k, x_hat_k_k1);
    Eigen::Matrix<double, dim_state, dim_state> F_k1;
    system_->jacobx_state_fcn(x_hat_k_k_, u_k, F_k1);
    Eigen::Matrix<double, dim_state, dim_state> P_k_k1 = F_k1 * P_k1_k1 * F_k1.transpose() + W_k1;

    // UPDATE
    Eigen::Matrix<double, dim_output, 1> y_hat_k_k1;
    system_->output_fcn(x_hat_k_k1, u_k, y_hat_k_k1);
    Eigen::Matrix<double, dim_output, 1> y_tilde_k = y_k - y_hat_k_k1;

    Eigen::Matrix<double, dim_output, dim_state> H_k;
    system_->jacobx_output_fcn(x_hat_k_k1, u_k, H_k);
    Eigen::Matrix<double, dim_output, dim_output> S_k = H_k * P_k_k1 * H_k.transpose() + V_k;
    Eigen::Matrix<double, dim_state, dim_output> K_k = P_k_k1 * H_k.transpose() * S_k.inverse();
    x_hat_k_k_ = x_hat_k_k1 + K_k * y_tilde_k;
    P_ = (Identity_x_ - K_k * H_k) * P_k_k1;
    system_->output_fcn(x_hat_k_k_, u_k, y_hat_k_);
  }

  Eigen::Matrix<double, dim_state, 1> get_state() const
  {
    return x_hat_k_k_;
  }
  Eigen::Matrix<double, dim_output, 1> get_output() const
  {
    return y_hat_k_;
  }

  void set_state(const Eigen::Ref<const Eigen::Matrix<double, dim_state, 1>>& x)
  {
    x_hat_k_k_ = x;
  }

  void reset()
  {
    P_ = W_;
    x_hat_k_k_.setZero();
    y_hat_k_.setZero();
  }

  void display()
  {
    std::cout << "\n" << std::endl;
    std::cout << "Extended Kalman Filter\n" << std::endl;
    std::cout << "W: \n" << W_ << std::endl;
    std::cout << "V: \n" << V_ << std::endl;
    system_->display();
  }

private:
  typename StateSpaceInterface::SharedPtr system_;
  Eigen::Matrix<double, dim_state, dim_state> W_;
  Eigen::Matrix<double, dim_output, dim_output> V_;
  Eigen::Matrix<double, dim_state, dim_state> P_;
  Eigen::Matrix<double, dim_state, dim_state> Identity_x_;
  Eigen::Matrix<double, dim_state, 1> x_hat_k_k_;
  Eigen::Matrix<double, dim_output, 1> y_hat_k_;
};

}  // namespace uclv::systems