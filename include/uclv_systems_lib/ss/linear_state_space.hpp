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

template <int dim1_state, int dim_input, int dim_output, int dim2_state = 1>
class LinearStateSpace
  : public StateSpaceInterface<dim1_state, dim_input, dim_output, dim2_state, dim2_state, dim2_state>
{
public:
  typedef std::shared_ptr<LinearStateSpace> SharedPtr;
  typedef std::shared_ptr<const LinearStateSpace> ConstSharedPtr;
  typedef std::weak_ptr<LinearStateSpace> WeakPtr;
  typedef std::weak_ptr<const LinearStateSpace> ConstWeakPtr;
  typedef std::unique_ptr<LinearStateSpace> UniquePtr;

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

  inline virtual const Eigen::Matrix<double, dim1_state, dim2_state>& get_state() const
  {
    return x_;
  }

  inline virtual const Eigen::Matrix<double, dim_output, dim2_state>& get_output() const
  {
    return y_;
  }

  /*==============================================*/

  /*=============SETTER===========================*/

  inline virtual void set_state(const Eigen::Ref<const Eigen::Matrix<double, dim1_state, dim2_state>>& x)
  {
    x_ = x;
    y_ = C * x_;
  }

  /*==============================================*/

  /*=============RUNNER===========================*/

  inline virtual const Eigen::Matrix<double, dim_output, dim2_state>&
  step(const Eigen::Ref<const Eigen::Matrix<double, dim_input, dim2_state>>& u_k)
  {
    x_ = A * x_ + B * u_k;
    y_ = C * x_ + D * u_k;
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
    std::cout << "Dim (input-state-output): " << dim_input << "x" << dim2_state << " - " << dim1_state << "x"
              << dim2_state << " - " << dim_output << "x" << dim2_state << "\n";
    std::cout << "A:\n" << A << "\n";
    std::cout << "B:\n" << B << "\n";
    std::cout << "C:\n" << C << "\n";
    std::cout << "D:\n" << D << "\n";
    std::cout << "current state:\n" << x_ << "\n";
  }

  /*==============================================*/

protected:
  Eigen::Matrix<double, dim1_state, dim2_state> x_;
  Eigen::Matrix<double, dim_output, dim2_state> y_;

public:
  Eigen::Matrix<double, dim1_state, dim1_state> A;
  Eigen::Matrix<double, dim1_state, dim_input> B;
  Eigen::Matrix<double, dim_output, dim1_state> C;
  Eigen::Matrix<double, dim_output, dim_input> D;
};

}  // namespace uclv::systems
