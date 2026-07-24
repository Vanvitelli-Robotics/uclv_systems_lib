/*
    Linear State Space System interface Class Discrete Time System

    Copyright 2024 Università della Campania Luigi Vanvitelli

    Author: Marco Costanzo <marco.costanzo@unicampania.it>

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

#include "state_space_interface.hpp"

/*! \file linear_state_space.hpp
    \brief This class represents a generic Discrete Time State Space System.
*/

namespace uclv::systems
{

template <typename Scalar_t, int dim_state, int dim_input, int dim_output, int num_col = 1>
class LinearStateSpace : public StateSpaceInterface<Scalar_t, dim_state, dim_input, dim_output, num_col, num_col, num_col>
{

  // TODO: Add support for num_col > 1 (for now, only num_col == 1 is supported)
  static_assert(num_col == 1, "LinearStateSpace supports only num_col == 1 (for now)");

public:

  using SharedPtr = std::shared_ptr<LinearStateSpace>;
  using ConstSharedPtr = std::shared_ptr<const LinearStateSpace>;
  using WeakPtr = std::weak_ptr<LinearStateSpace>;
  using ConstWeakPtr = std::weak_ptr<const LinearStateSpace>;
  using UniquePtr = std::unique_ptr<LinearStateSpace>;

  using StateSpaceInterface_t = StateSpaceInterface<Scalar_t, dim_state, dim_input, dim_output, num_col, num_col, num_col>;
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

  using DynamicMatrix_t = Eigen::Matrix<Scalar_t, dim_state, dim_state>;
  using DynamicMatrixRef_t = Eigen::Ref<DynamicMatrix_t>;
  using DynamicMatrixConstRef_t = Eigen::Ref<const DynamicMatrix_t>;

  using InputMatrix_t = Eigen::Matrix<Scalar_t, dim_state, dim_input>;
  using InputMatrixRef_t = Eigen::Ref<InputMatrix_t>;
  using InputMatrixConstRef_t = Eigen::Ref<const InputMatrix_t>;

  using OutputMatrix_t = Eigen::Matrix<Scalar_t, dim_output, dim_state>;
  using OutputMatrixRef_t = Eigen::Ref<OutputMatrix_t>;
  using OutputMatrixConstRef_t = Eigen::Ref<const OutputMatrix_t>;

  using FeedthroughMatrix_t = Eigen::Matrix<Scalar_t, dim_output, dim_input>;
  using FeedthroughMatrixRef_t = Eigen::Ref<FeedthroughMatrix_t>;
  using FeedthroughMatrixConstRef_t = Eigen::Ref<const FeedthroughMatrix_t>;

  /*===============CONSTRUCTORS===================*/

  LinearStateSpace() = default;

  //! Copy Constructor
  LinearStateSpace(const LinearStateSpace& sys) = default;

  virtual ~LinearStateSpace() = default;

  //! Clone the object
  virtual LinearStateSpace* clone() const
  {
    return new LinearStateSpace(*this);
  }

  /*==============================================*/

  /*=============GETTER===========================*/

  inline virtual const State_t& get_state() const
  {
    return x_;
  }

  inline virtual const Output_t& get_output() const
  {
    return y_;
  }

  /*==============================================*/

  /*=============SETTER===========================*/

  inline virtual void set_state(const StateConstRef_t& x)
  {
    x_ = x;
    y_ = C * x_;
  }

  /*==============================================*/

  /*=============RUNNER===========================*/
  inline virtual void state_fcn(const StateConstRef_t& x,
                                const InputConstRef_t& u_k,
                                StateRef_t out) const
  {
    out = A * x + B * u_k;
  }
  inline virtual void output_fcn(const StateConstRef_t& x,
                                 const InputConstRef_t& u_k,
                                 OutputRef_t out) const
  {
    out = C * x + D * u_k;
  }

  inline virtual void jacobx_state_fcn(const StateConstRef_t& x,
                                       const InputConstRef_t& u_k,
                                       JacobianStateStateRef_t out) const
  {
    (void)x;
    (void)u_k;
    if (num_col != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for num_col != 1");
    }
    out = A;
  }
  inline virtual void jacobu_state_fcn(const StateConstRef_t& x,
                                       const InputConstRef_t& u_k,
                                       JacobianStateInputRef_t out) const
  {
    (void)x;
    (void)u_k;
    if (num_col != 1)
    {
      throw std::runtime_error("The Jacobian of the state function is not defined for num_col != 1");
    }
    out = B;
  }

  inline virtual void jacobx_output_fcn(const StateConstRef_t& x,
                                        const InputConstRef_t& u_k,
                                        JacobianOutputStateRef_t out) const
  {
    (void)x;
    (void)u_k;
    if (num_col != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for num_col != 1");
    }
    out = C;
  }

  inline virtual void jacobu_output_fcn(const StateConstRef_t& x,
                                        const InputConstRef_t& u_k,
                                        JacobianOutputInputRef_t out) const
  {
    (void)x;
    (void)u_k;
    if (num_col != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for num_col != 1");
    }
    out = D;
  }

  inline virtual const Output_t& step(const InputConstRef_t& u_k)
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
    std::cout << "Linear State Space System\n";
    std::cout << "Dim (input-state-output): " << dim_input << "x" << num_col << " - " << dim_state << "x" << num_col
              << " - " << dim_output << "x" << num_col << "\n";
    std::cout << "A:\n" << A << "\n";
    std::cout << "B:\n" << B << "\n";
    std::cout << "C:\n" << C << "\n";
    std::cout << "D:\n" << D << "\n";
    std::cout << "current state:\n" << x_ << "\n";
  }

  /*==============================================*/

protected:
  State_t x_;
  Output_t y_;

public:
  DynamicMatrix_t A;
  InputMatrix_t B;
  OutputMatrix_t C;
  FeedthroughMatrix_t D;
};

}  // namespace uclv::systems
