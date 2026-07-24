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

  template <typename Scalar_t, int dim1_state, int dim1_input, int dim1_output, int dim2_state = 1, int dim2_input = 1, int dim2_output = 1, typename size_t = std::size_t>
  class StateSpaceInterface : public SystemInterface<Scalar_t, dim1_input, dim1_output, dim2_input, dim2_output, size_t>
  {
  public:

    using SharedPtr = std::shared_ptr<StateSpaceInterface>;
    using ConstSharedPtr = std::shared_ptr<const StateSpaceInterface>;
    using WeakPtr = std::weak_ptr<StateSpaceInterface>;
    using ConstWeakPtr = std::weak_ptr<const StateSpaceInterface>;
    using UniquePtr = std::unique_ptr<StateSpaceInterface>;

    using SystemInterface_t = SystemInterface<Scalar_t, dim1_input, dim1_output, dim2_input, dim2_output, size_t>;
    using Input_t = typename SystemInterface_t::Input_t;
    using InputRef_t = typename SystemInterface_t::InputRef_t;
    using InputConstRef_t = typename SystemInterface_t::InputConstRef_t;
    using Output_t = typename SystemInterface_t::Output_t;
    using OutputRef_t = typename SystemInterface_t::OutputRef_t;
    using OutputConstRef_t = typename SystemInterface_t::OutputConstRef_t;

    using State_t = Eigen::Matrix<Scalar_t, dim1_state, dim2_state>;
    using StateRef_t = Eigen::Ref<State_t>;
    using StateConstRef_t = Eigen::Ref<const State_t>;

    using JacobianStateState_t = Eigen::Matrix<Scalar_t, dim1_state * dim2_state, dim1_state * dim2_state>;
    using JacobianStateStateRef_t = Eigen::Ref<JacobianStateState_t>;
    using JacobianStateStateConstRef_t = Eigen::Ref<const JacobianStateState_t>;

    using JacobianStateInput_t = Eigen::Matrix<Scalar_t, dim1_state * dim2_state, dim1_input * dim2_input>;
    using JacobianStateInputRef_t = Eigen::Ref<JacobianStateInput_t>;
    using JacobianStateInputConstRef_t = Eigen::Ref<const JacobianStateInput_t>;

    using JacobianOutputState_t = Eigen::Matrix<Scalar_t, dim1_output * dim2_output, dim1_state * dim2_state>;
    using JacobianOutputStateRef_t = Eigen::Ref<JacobianOutputState_t>;
    using JacobianOutputStateConstRef_t = Eigen::Ref<const JacobianOutputState_t>;

    using JacobianOutputInput_t = Eigen::Matrix<Scalar_t, dim1_output * dim2_output, dim1_input * dim2_input>;
    using JacobianOutputInputRef_t = Eigen::Ref<JacobianOutputInput_t>;
    using JacobianOutputInputConstRef_t = Eigen::Ref<const JacobianOutputInput_t>;

  protected:
  public:
    /*===============CONSTRUCTORS===================*/

    StateSpaceInterface() = default;

    //! Copy Constructor
    StateSpaceInterface(const StateSpaceInterface &sys) = default;

    virtual ~StateSpaceInterface() = default;

    //! Clone the object
    virtual StateSpaceInterface *clone() const = 0;

    /*==============================================*/

    /*=============GETTER===========================*/

    inline virtual const State_t &get_state() const = 0;

    inline virtual const Output_t &get_output() const = 0;

    /*==============================================*/

    /*=============SETTER===========================*/

    inline virtual void set_state(const StateConstRef_t &x) = 0;

    /*==============================================*/

    /*=============RUNNER===========================*/

    //! State function
    inline virtual void state_fcn(const StateConstRef_t &x,
                                  const InputConstRef_t &u_k,
                                  StateRef_t out) const = 0;

    //! Output function
    inline virtual void output_fcn(const StateConstRef_t &x,
                                   const InputConstRef_t &u_k,
                                   OutputRef_t out) const = 0;

    //! Jacobian of the state function with respect to the state
    inline virtual void jacobx_state_fcn(const StateConstRef_t &x,
                                         const InputConstRef_t &u_k,
                                         JacobianStateStateRef_t out) const
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
    inline virtual void jacobu_state_fcn(const StateConstRef_t &x,
                                         const InputConstRef_t &u_k,
                                         JacobianStateInputRef_t out) const
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
    inline virtual void jacobx_output_fcn(const StateConstRef_t &x,
                                          const InputConstRef_t &u_k,
                                          JacobianOutputStateRef_t out) const
    {
      (void)x;
      (void)u_k;
      (void)out;
      if (dim2_output != 1)
      {
        throw std::runtime_error("The Jacobian of the output function is not defined for dim2_output != 1");
      }
    }

    //! Jacobian of the output function with respect to the input
    inline virtual void jacobu_output_fcn(const StateConstRef_t &x,
                                          const InputConstRef_t &u_k,
                                          JacobianOutputInputRef_t out) const
    {
      (void)x;
      (void)u_k;
      (void)out;
      if (dim2_output != 1)
      {
        throw std::runtime_error("The Jacobian of the output function is not defined for dim2_output != 1");
      }
    }

    inline virtual const Output_t &
    step(const InputConstRef_t &u_k) override = 0;

    /*==============================================*/

    /*=============VARIE===========================*/
    inline virtual void reset() = 0;

    inline virtual void get_resetted_state(StateRef_t x) const 
    {
      x.setZero();
    }

    inline virtual void get_resetted_output(OutputRef_t y) const
    {
      y.setZero();
    }

    virtual size_t get_size_state() const
    {
      return dim1_state;
    }

    virtual size_t get_size1_state() const
    {
      return dim1_state;
    }

    virtual size_t get_size2_state() const
    {
      return dim2_state;
    }

    virtual void display() const = 0;

    /*==============================================*/
  };

} // namespace uclv::systems
