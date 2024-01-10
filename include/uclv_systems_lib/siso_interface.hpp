/*
    SISO interface Class Discrete Time Transfer Function SISO (Linear Filter)

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

/*! \file siso.hpp
    \brief This class represents a generic Discrete Time SISO System.
*/

namespace uclv::systems
{

class SISOInterface
{
public:
  typedef std::shared_ptr<SISOInterface> SharedPtr;
  typedef std::shared_ptr<const SISOInterface> ConstSharedPtr;
  typedef std::weak_ptr<SISOInterface> WeakPtr;
  typedef std::weak_ptr<const SISOInterface> ConstWeakPtr;
  typedef std::unique_ptr<SISOInterface> UniquePtr;

protected:
public:
  /*===============CONSTRUCTORS===================*/

  SISOInterface() = default;

  //! Copy Constructor
  SISOInterface(const SISOInterface& tf) = default;

  virtual ~SISOInterface() = default;

  //! Clone the object
  virtual SISOInterface* clone() const = 0;

  /*==============================================*/

  /*=============GETTER===========================*/

  /*==============================================*/

  /*=============SETTER===========================*/

  /*==============================================*/

  /*=============RUNNER===========================*/

  inline virtual double step(double u_k) = 0;

  /*==============================================*/

  /*=============VARIE===========================*/
  inline virtual void reset() = 0;

  virtual unsigned int getSizeInput() const
  {
    return 1;
  }

  virtual unsigned int getSizeOutput() const
  {
    return 1;
  }

  //! Get the last output
  virtual double getLastOutput() const = 0;

  virtual void display() const = 0;

  /*==============================================*/
};

}  // namespace uclv::systems
