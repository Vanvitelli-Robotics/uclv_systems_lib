/*
    System interface Class Discrete Time System

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

// TODO for the future:
//  In C++ moderno e in Eigen si preferisce il CRTP (Curiously Recurring Template Pattern) per fare polimorfismo a tempo di compilazione senza metodi virtuali.
// Use the CRTP pattern to avoid virtual functions and improve performance. This will require refactoring the current design to use templates and static polymorphism.

#pragma once

#include <memory>
#include <iostream>
#include <Eigen/Dense>

/*! \file system_interface.hpp
    \brief This class represents a generic Discrete Time System.
*/

namespace uclv::systems
{

  template <typename Scalar_t, int dim1_input, int dim1_output, int dim2_input = 1, int dim2_output = 1, typename size_t = std::size_t>
  class SystemInterface
  {
  public:
    using SharedPtr = std::shared_ptr<SystemInterface>;
    using ConstSharedPtr = std::shared_ptr<const SystemInterface>;
    using WeakPtr = std::weak_ptr<SystemInterface>;
    using ConstWeakPtr = std::weak_ptr<const SystemInterface>;
    using UniquePtr = std::unique_ptr<SystemInterface>;

    using Input_t = Eigen::Matrix<Scalar_t, dim1_input, dim2_input>;
    using InputRef_t = Eigen::Ref<Input_t>;
    using InputConstRef_t = Eigen::Ref<const Input_t>;
    using Output_t = Eigen::Matrix<Scalar_t, dim1_output, dim2_output>;
    using OutputRef_t = Eigen::Ref<Output_t>;
    using OutputConstRef_t = Eigen::Ref<const Output_t>;

  protected:
  public:
    /*===============CONSTRUCTORS===================*/

    SystemInterface() = default;

    //! Copy Constructor
    SystemInterface(const SystemInterface &sys) = default;

    virtual ~SystemInterface() = default;

    //! Clone the object
    virtual SystemInterface *clone() const = 0;

    /*==============================================*/

    /*=============GETTER===========================*/

    inline virtual const Output_t &get_output() const = 0;

    /*==============================================*/

    /*=============SETTER===========================*/

    /*==============================================*/

    /*=============RUNNER===========================*/

    inline virtual const Output_t &
    step(const InputConstRef_t &u_k) = 0;

    /*==============================================*/

    /*=============VARIE===========================*/
    inline virtual void reset() = 0;

    virtual size_t get_size_input() const
    {
      return dim1_input;
    }

    virtual size_t get_size1_input() const
    {
      return dim1_input;
    }

    virtual size_t get_size2_input() const
    {
      return dim2_input;
    }

    virtual size_t get_size_output() const
    {
      return dim1_output;
    }

    virtual size_t get_size1_output() const
    {
      return dim1_output;
    }

    virtual size_t get_size2_output() const
    {
      return dim2_output;
    }

    virtual void display() const = 0;

    /*==============================================*/
  };

} // namespace uclv::systems
