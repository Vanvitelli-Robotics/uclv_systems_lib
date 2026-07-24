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

#include "discretizator_interface.hpp"
#include "../continuous_time/continuous_time_state_space_interface.hpp"

/*! \file forward_euler.hpp
    \brief This class represents the Forward Euler Discretizator.
*/

namespace uclv::systems
{

template <typename Scalar_t, int dim1_state, int dim1_input, int dim1_output, int dim2_state = 1, int dim2_input = 1, int dim2_output = 1>
class ForwardEuler
  : public DiscretizatorInterface<Scalar_t, dim1_state, dim1_input, dim1_output, dim2_state, dim2_input, dim2_output>
{
public:
  typedef std::shared_ptr<ForwardEuler> SharedPtr;
  typedef std::shared_ptr<const ForwardEuler> ConstSharedPtr;
  typedef std::weak_ptr<ForwardEuler> WeakPtr;
  typedef std::weak_ptr<const ForwardEuler> ConstWeakPtr;
  typedef std::unique_ptr<ForwardEuler> UniquePtr;

  typedef ::uclv::systems::ContinuousTimeStateSpaceInterface<Scalar_t, dim1_state, dim1_input, dim1_output, dim2_state,
                                                             dim2_input, dim2_output>
      ContinuousTimeStateSpaceInterface;

  /*===============CONSTRUCTORS===================*/

  ForwardEuler() = default;

  ForwardEuler(typename ContinuousTimeStateSpaceInterface::SharedPtr system_ptr, Scalar_t sample_time)
    : sys_(system_ptr), sample_time_(sample_time)
  {
    x_.resizeLike(sys_->get_state());
    y_.resizeLike(sys_->get_output());
  }

  ForwardEuler(typename ContinuousTimeStateSpaceInterface::SharedPtr system_ptr, Scalar_t sample_time, const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, 1>>& x)
    : sys_(system_ptr), sample_time_(sample_time)
  {
    x_.resizeLike(sys_->get_state());
    x_ = x;
    y_.resizeLike(sys_->get_output());
    y_ = sys_->get_output();
  }

  //! Copy Constructor
  ForwardEuler(const ForwardEuler& sys) = default;

  virtual ~ForwardEuler() = default;

  //! Clone the object
  virtual ForwardEuler* clone() const
  {
    return new ForwardEuler(*this);
  }

  /*==============================================*/

  /*=============GETTER===========================*/

  inline virtual const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>& get_state() const
  {
    return x_;
  }

  inline virtual const Eigen::Matrix<Scalar_t, dim1_output, dim2_output>& get_output() const
  {
    return y_;
  }

  inline virtual Scalar_t get_sample_time() const
  {
    return sample_time_;
  }

  /*==============================================*/

  /*=============SETTER===========================*/

  inline virtual void set_state(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& x)
  {
    x_ = x;
  }

  /*==============================================*/

  /*=============RUNNER===========================*/
  inline virtual void state_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& x,
                                const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_input, dim2_input>>& u_k,
                                Eigen::Ref<Eigen::Matrix<Scalar_t, dim1_state, dim2_state>> out) const
  {
    Eigen::Matrix<Scalar_t, dim1_state, dim2_state> x_k1;
    sys_->state_fcn(x, u_k, x_k1);
    out = x + sample_time_ * x_k1;
  }
  inline virtual void output_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& x,
                                 const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_input, dim2_input>>& u_k,
                                 Eigen::Ref<Eigen::Matrix<Scalar_t, dim1_output, dim2_output>> out) const
  {
    sys_->output_fcn(x, u_k, out);
  }

  inline virtual void jacobx_state_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& x,
                                       const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_input, dim2_input>>& u_k,
                                       Eigen::Ref<Eigen::Matrix<Scalar_t, dim1_state, dim1_state>> out) const
  {
    (void)x;
    (void)u_k;
    if (dim2_state != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for dim2_state != 1");
    }
    sys_->jacobx_state_fcn(x, u_k, out);
    out = Eigen::Matrix<Scalar_t, dim1_state, dim1_state>::Identity() + sample_time_ * out;
  }
  inline virtual void jacobu_state_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& x,
                                       const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_input, dim2_input>>& u_k,
                                       Eigen::Ref<Eigen::Matrix<Scalar_t, dim1_state, dim1_input>> out) const
  {
    (void)x;
    (void)u_k;
    if (dim2_state != 1)
    {
      throw std::runtime_error("The Jacobian of the state function is not defined for dim2_state != 1");
    }
    sys_->jacobu_state_fcn(x, u_k, out);
  }

  inline virtual void jacobx_output_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& x,
                                        const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_input, dim2_input>>& u_k,
                                        Eigen::Ref<Eigen::Matrix<Scalar_t, dim1_output, dim1_state>> out) const
  {
    (void)x;
    (void)u_k;
    if (dim2_output != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for dim2_output != 1");
    }
    sys_->jacobx_output_fcn(x, u_k, out);
  }

  inline virtual void jacobu_output_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& x,
                                        const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_input, dim2_input>>& u_k,
                                        Eigen::Ref<Eigen::Matrix<Scalar_t, dim1_output, dim1_input>> out) const
  {
    (void)x;
    (void)u_k;
    if (dim2_output != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for dim2_output != 1");
    }
    sys_->jacobu_output_fcn(x, u_k, out);
  }

  inline virtual const Eigen::Matrix<Scalar_t, dim1_output, dim2_output>&
  step(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_input, dim2_input>>& u_k)
  {
    state_fcn(x_, u_k, x_);
    output_fcn(x_, u_k, y_);
    return y_;
  }

  /*==============================================*/

  /*=============VARIE===========================*/
  inline virtual void reset()
  {
    x_.setZero();
    y_.setZero();
  }

  virtual void display() const
  {
    std::cout << "Forward Euler Discretizator" << std::endl;
    std::cout << "Discretizator State: " << std::endl;
    std::cout << x_ << std::endl;
    std::cout << "Sample Time: " << sample_time_ << std::endl;
    sys_->display();
  }

  /*==============================================*/

protected:
  typename ContinuousTimeStateSpaceInterface::SharedPtr sys_;
  Scalar_t sample_time_;
  Eigen::Matrix<Scalar_t, dim1_state, dim2_state> x_;
  Eigen::Matrix<Scalar_t, dim1_output, dim2_output> y_;
};

}  // namespace uclv::systems
