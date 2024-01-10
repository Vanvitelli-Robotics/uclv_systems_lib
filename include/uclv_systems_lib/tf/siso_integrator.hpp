/*
    TF_INTEGTATOR Class

    Integrator transfer function using the trapez method

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

/*! \file TF_INTEGRATOR.h
    \brief This class represents an integrator as Discrete Time Transfer Function System.
*/

#include "siso.hpp"

namespace uclv::systems::tf::siso
{
//!  TF_INTEGRATOR class: represents an integrator as Discrete Time Transfer Function System.
/*!
    This class is an integrator as Discrete Time Transfer Function System.

    It is discretized using the tustin method.

    \sa Linear_System_Interface, TF_FIRST_ORDER_FILTER, TF_SISO, TF_MIMO
*/
class Integrator : public SISO
{
private:
  Integrator();  // NO DEFAULT CONSTRUCTOR

protected:
  double Ts_;  //!< Sampling time

  //! Integrator gain
  double gain_;

public:
  typedef std::shared_ptr<Integrator> SharedPtr;
  typedef std::shared_ptr<const Integrator> ConstSharedPtr;
  typedef std::weak_ptr<Integrator> WeakPtr;
  typedef std::weak_ptr<const Integrator> ConstWeakPtr;
  typedef std::unique_ptr<Integrator> UniquePtr;

  /*===============CONSTRUCTORS===================*/

  //! Constructor
  /*!
      \param Ts sampling time
      \param gain gain of the integrator, default = 1
  */
  Integrator(double Ts, double gain = 1.0) : SISO(compute_num_coeff(Ts), compute_den_coeff()), Ts_(Ts), gain_(gain)
  {
  }

  //! Destructor
  virtual ~Integrator() override = default;

  virtual Integrator* clone() const override
  {
    return new Integrator(*this);
  }

  //! Copy Constructor
  Integrator(const Integrator& tf) = default;
  /*==============================================*/

  /*===============STATIC FUNCTIONS FOR SIMPLE CONSTRUCTOR WRITING=======*/

  //! INTERNAL, numerator coefficients
  static Eigen::Vector2d compute_num_coeff(double Ts)
  {
    return (Ts / 2.0) * Eigen::Vector2d::Ones();
  }

  //! INTERNAL, denominator coefficients
  static Eigen::Vector2d compute_den_coeff()
  {
    return Eigen::Vector2d(1.0, -1.0);
  }

  /*==============================================*/

  /*=============GETTER===========================*/
  /*==============================================*/

  /*=============SETTER===========================*/
  inline virtual void setTs(double Ts)
  {
    Ts_ = Ts;
    b_vec_ = compute_num_coeff(Ts_);
  }

  //! Change the state so that the output is "output"
  /*!
      \param output output after the state change
  */
  inline virtual void setOutput(double output)
  {
    u_vec_[0] = 0.0;
    u_vec_[1] = 0.0;
    y_vec_[0] = output;
    y_k_ = output;
  }
  /*==============================================*/

  /*=============RUNNER===========================*/
  inline virtual double step(double uk) override
  {
    return SISO::step(gain_ * uk);
  }
  /*==============================================*/

  /*=============VARIE===========================*/

  virtual void display_tf() const override
  {
    std::cout << "TF_INTEGRATOR:" << std::endl << "   Ts: " << Ts_ << std::endl << "   gain: " << gain_ << std::endl;
  }

  /*==============================================*/
};

/*=============STATIC FUNS===========================*/
/*==============================================*/

}  // namespace uclv::systems::tf::siso
