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

#pragma once

#include <memory>
#include <iostream>
#include <Eigen/Dense>

/*! \file system_interface.hpp
    \brief This class represents a generic Discrete Time System.
*/

namespace uclv::systems
{

template <int dim1_input, int dim1_output, int dim2_input = 1, int dim2_output = 1>
class SystemInterface
{
public:
  typedef std::shared_ptr<SystemInterface> SharedPtr;
  typedef std::shared_ptr<const SystemInterface> ConstSharedPtr;
  typedef std::weak_ptr<SystemInterface> WeakPtr;
  typedef std::weak_ptr<const SystemInterface> ConstWeakPtr;
  typedef std::unique_ptr<SystemInterface> UniquePtr;

protected:
public:
  /*===============CONSTRUCTORS===================*/

  SystemInterface() = default;

  //! Copy Constructor
  SystemInterface(const SystemInterface& sys) = default;

  virtual ~SystemInterface() = default;

  //! Clone the object
  virtual SystemInterface* clone() const = 0;

  /*==============================================*/

  /*=============GETTER===========================*/

  inline virtual const Eigen::Matrix<double, dim1_output, dim2_output>& get_output() const = 0;

  /*==============================================*/

  /*=============SETTER===========================*/

  /*==============================================*/

  /*=============RUNNER===========================*/

  inline virtual const Eigen::Matrix<double, dim1_output, dim2_output>&
  step(const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k) = 0;

  /*==============================================*/

  /*=============VARIE===========================*/
  inline virtual void reset() = 0;

  virtual unsigned int get_size_input() const
  {
    return dim1_input;
  }

  virtual unsigned int get_size1_input() const
  {
    return dim1_input;
  }

  virtual unsigned int get_size2_input() const
  {
    return dim2_input;
  }

  virtual unsigned int get_size_output() const
  {
    return dim1_output;
  }

  virtual unsigned int get_size1_output() const
  {
    return dim1_output;
  }

  virtual unsigned int get_size2_output() const
  {
    return dim2_output;
  }

  virtual void display() const = 0;

  /*==============================================*/
};

}  // namespace uclv::systems
