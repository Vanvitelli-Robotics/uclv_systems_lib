/*
    State Space System interface Class Discrete Time System

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

#include "../system_interface.hpp"

/*! \file state_space_interface.hpp
    \brief This class represents a generic Discrete Time State Space System.
*/

namespace uclv::systems
{

template <int dim1_state, int dim1_input, int dim1_output, int dim2_state = 1, int dim2_input = 1, int dim2_output = 1>
class StateSpaceInterface : public SystemInterface<dim1_input, dim1_output, dim2_input, dim2_output>
{
public:
  typedef std::shared_ptr<StateSpaceInterface> SharedPtr;
  typedef std::shared_ptr<const StateSpaceInterface> ConstSharedPtr;
  typedef std::weak_ptr<StateSpaceInterface> WeakPtr;
  typedef std::weak_ptr<const StateSpaceInterface> ConstWeakPtr;
  typedef std::unique_ptr<StateSpaceInterface> UniquePtr;

protected:
public:
  /*===============CONSTRUCTORS===================*/

  StateSpaceInterface() = default;

  //! Copy Constructor
  StateSpaceInterface(const StateSpaceInterface& sys) = default;

  virtual ~StateSpaceInterface() = default;

  //! Clone the object
  virtual StateSpaceInterface* clone() const = 0;

  /*==============================================*/

  /*=============GETTER===========================*/

  inline virtual const Eigen::Matrix<double, dim1_state, dim2_state>& get_state() const = 0;

  inline virtual const Eigen::Matrix<double, dim1_output, dim2_output>& get_output() const = 0;

  /*==============================================*/

  /*=============SETTER===========================*/

  inline virtual void set_state(const Eigen::Ref<const Eigen::Matrix<double, dim1_state, dim2_state>>& x) = 0;

  /*==============================================*/

  /*=============RUNNER===========================*/

  //! State function
  inline virtual void state_fcn(const Eigen::Ref<const Eigen::Matrix<double, dim1_state, dim2_state>>& x,
                                const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k,
                                Eigen::Ref<Eigen::Matrix<double, dim1_state, dim2_state>> out) const = 0;

  //! Output function
  inline virtual void output_fcn(const Eigen::Ref<const Eigen::Matrix<double, dim1_state, dim2_state>>& x,
                                 const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k,
                                 Eigen::Ref<Eigen::Matrix<double, dim1_output, dim2_output>> out) const = 0;

  //! Jacobian of the state function with respect to the state
  inline virtual void jacobx_state_fcn(const Eigen::Ref<const Eigen::Matrix<double, dim1_state, dim2_state>>& x,
                                       const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k,
                                       Eigen::Ref<Eigen::Matrix<double, dim1_state, dim1_state>> out) const
  {
    (void)x;
    (void)u_k;
    (void)out;
    if (dim2_state != 1)
    {
      throw std::runtime_error("The Jacobian of the state function is not defined for dim2_state != 1");
    }
  }

  //! Jacobian of the state function with respect to the input
  inline virtual void jacobu_state_fcn(const Eigen::Ref<const Eigen::Matrix<double, dim1_state, dim2_state>>& x,
                                       const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k,
                                       Eigen::Ref<Eigen::Matrix<double, dim1_state, dim1_input>> out) const
  {
    (void)x;
    (void)u_k;
    (void)out;
    if (dim2_state != 1)
    {
      throw std::runtime_error("The Jacobian of the state function is not defined for dim2_state != 1");
    }
  }

  //! Jacobian of the output function with respect to the state
  inline virtual void jacobx_output_fcn(const Eigen::Ref<const Eigen::Matrix<double, dim1_state, dim2_state>>& x,
                                        const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k,
                                        Eigen::Ref<Eigen::Matrix<double, dim1_output, dim1_state>> out) const
  {
    (void)x;
    (void)u_k;
    (void)out;
    if (dim2_output != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for dim2_output != 1");
    }
  }

  // Jacobian of the output function with respect to the input
  inline virtual void jacobu_output_fcn(const Eigen::Ref<const Eigen::Matrix<double, dim1_state, dim2_state>>& x,
                                        const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k,
                                        Eigen::Ref<Eigen::Matrix<double, dim1_output, dim1_input>> out) const
  {
    (void)x;
    (void)u_k;
    (void)out;
    if (dim2_output != 1)
    {
      throw std::runtime_error("The Jacobian of the output function is not defined for dim2_output != 1");
    }
  }

  inline virtual const Eigen::Matrix<double, dim1_output, dim2_output>&
  step(const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k) = 0;

  /*==============================================*/

  /*=============VARIE===========================*/
  inline virtual void reset() = 0;

  virtual unsigned int get_size_state() const
  {
    return dim1_state;
  }

  virtual unsigned int get_size1_state() const
  {
    return dim1_state;
  }

  virtual unsigned int get_size2_state() const
  {
    return dim2_state;
  }

  virtual void display() const = 0;

  /*==============================================*/
};

}  // namespace uclv::systems
