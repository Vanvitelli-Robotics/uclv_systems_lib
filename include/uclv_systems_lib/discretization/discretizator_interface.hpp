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

#include "../ss/state_space_interface.hpp"

/*! \file discretizator_interface.hpp
    \brief This class represents a generic Discretizator.
*/

namespace uclv::systems
{

template <typename Scalar_t, int dim1_state, int dim1_input, int dim1_output, int dim2_state = 1, int dim2_input = 1, int dim2_output = 1>
class DiscretizatorInterface
  : public StateSpaceInterface<Scalar_t, dim1_state, dim1_input, dim1_output, dim2_state, dim2_input, dim2_output>
{
public:
  typedef std::shared_ptr<DiscretizatorInterface> SharedPtr;
  typedef std::shared_ptr<const DiscretizatorInterface> ConstSharedPtr;
  typedef std::weak_ptr<DiscretizatorInterface> WeakPtr;
  typedef std::weak_ptr<const DiscretizatorInterface> ConstWeakPtr;
  typedef std::unique_ptr<DiscretizatorInterface> UniquePtr;

  /*===============CONSTRUCTORS===================*/

  DiscretizatorInterface() = default;

  //! Copy Constructor
  DiscretizatorInterface(const DiscretizatorInterface& sys) = default;

  virtual ~DiscretizatorInterface() = default;

  //! Clone the object
  virtual DiscretizatorInterface* clone() const = 0;

  /*==============================================*/

  /*=============GETTER===========================*/
  inline virtual double get_sample_time() const = 0;
};

}  // namespace uclv::systems
