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

#include "continuous_time_state_space_interface.hpp"

/*! \file linear_state_space.hpp
    \brief This class represents a generic Discrete Time State Space System.
*/

namespace uclv::systems
{

template <typename Scalar_t, int dim_state, int dim_input, int dim_output, int num_col = 1>
class ContinuousTimeLinearStateSpace
  : public ContinuousTimeStateSpaceInterface<Scalar_t, dim_state, dim_input, dim_output, num_col, num_col, num_col>
{
public:
  typedef std::shared_ptr<ContinuousTimeLinearStateSpace> SharedPtr;
  typedef std::shared_ptr<const ContinuousTimeLinearStateSpace> ConstSharedPtr;
  typedef std::weak_ptr<ContinuousTimeLinearStateSpace> WeakPtr;
  typedef std::weak_ptr<const ContinuousTimeLinearStateSpace> ConstWeakPtr;
  typedef std::unique_ptr<ContinuousTimeLinearStateSpace> UniquePtr;

  /*===============CONSTRUCTORS===================*/

  ContinuousTimeLinearStateSpace(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, dim_state>>& A,
                                 const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, dim_input>>& B,
                                 const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_output, dim_state>>& C,
                                 const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_output, dim_input>>& D,
                                 const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, num_col>>& x0)
    : A(A), B(B), C(C), D(D)
  {
    set_state(x0);
  }

  ContinuousTimeLinearStateSpace() = default;

  //! Copy Constructor
  ContinuousTimeLinearStateSpace(const ContinuousTimeLinearStateSpace& sys) = default;

  virtual ~ContinuousTimeLinearStateSpace() = default;

  //! Clone the object
  virtual ContinuousTimeLinearStateSpace* clone() const
  {
    return new ContinuousTimeLinearStateSpace(*this);
  }

  /*==============================================*/

  /*=============GETTER===========================*/

  inline virtual const Eigen::Matrix<Scalar_t, dim_state, num_col>& get_state() const
  {
    return x_;
  }

  inline virtual const Eigen::Matrix<Scalar_t, dim_output, num_col>& get_output() const
  {
    return y_;
  }

  /*==============================================*/

  /*=============SETTER===========================*/

  inline virtual void set_state(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, num_col>>& x)
  {
    x_ = x;
    y_ = C * x_;
  }

  /*==============================================*/

  /*=============RUNNER===========================*/
  inline virtual void state_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, num_col>>& x,
                                const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_input, num_col>>& u_k,
                                Eigen::Matrix<Scalar_t, dim_state, num_col>& out) const
  {
    out = A * x + B * u_k;
  }
  inline virtual void output_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, num_col>>& x,
                                 const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_input, num_col>>& u_k,
                                 Eigen::Matrix<Scalar_t, dim_output, num_col>& out) const
  {
    out = C * x + D * u_k;
  }

  inline virtual void jacobx_state_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, num_col>>& x,
                                       const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_input, num_col>>& u_k,
                                       Eigen::Matrix<Scalar_t, dim_state, dim_state>& out) const
  {
    (void)x;
    (void)u_k;
    if (num_col != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for num_col != 1");
    }
    out = A;
  }
  inline virtual void jacobu_state_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, num_col>>& x,
                                       const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_input, num_col>>& u_k,
                                       Eigen::Matrix<Scalar_t, dim_state, dim_input>& out) const
  {
    (void)x;
    (void)u_k;
    if (num_col != 1)
    {
      throw std::runtime_error("The Jacobian of the state function is not defined for num_col != 1");
    }
    out = B;
  }

  inline virtual void jacobx_output_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, num_col>>& x,
                                        const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_input, num_col>>& u_k,
                                        Eigen::Matrix<Scalar_t, dim_output, dim_state>& out) const
  {
    (void)x;
    (void)u_k;
    if (num_col != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for num_col != 1");
    }
    out = C;
  }

  inline virtual void jacobu_output_fcn(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_state, num_col>>& x,
                                        const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim_input, num_col>>& u_k,
                                        Eigen::Matrix<Scalar_t, dim_output, dim_input>& out) const
  {
    (void)x;
    (void)u_k;
    if (num_col != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for num_col != 1");
    }
    out = D;
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
    std::cout << "Continuous Time Linear State Space System\n";
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
  Eigen::Matrix<Scalar_t, dim_state, num_col> x_;
  Eigen::Matrix<Scalar_t, dim_output, num_col> y_;

public:
  Eigen::Matrix<Scalar_t, dim_state, dim_state> A;
  Eigen::Matrix<Scalar_t, dim_state, dim_input> B;
  Eigen::Matrix<Scalar_t, dim_output, dim_state> C;
  Eigen::Matrix<Scalar_t, dim_output, dim_input> D;
};

}  // namespace uclv::systems
