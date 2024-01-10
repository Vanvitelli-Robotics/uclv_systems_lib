/*
    PI CONTROLLER Class

    PI CONTROLLER transfer function using the trapez method

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

#include "../tf/siso.hpp"
#include "../tf/siso_integrator.hpp"

namespace uclv::systems::controllers
{
//!  PI CONTROLLER class: represents an integrator as Discrete Time Transfer Function System.
/*!
    This class is an integrator as Discrete Time Transfer Function System.

    It is discretized using the tustin method.

    \sa Linear_System_Interface, TF_FIRST_ORDER_FILTER, TF_SISO, TF_MIMO
*/
class PI : public SISOInterface
{
private:
  PI();  // NO DEFAULT CONSTRUCTOR

protected:
  double Ts_;  //!< Sampling time

  //! Integrator gain
  double Kp_, Ki_;

  double last_output_;

  uclv::systems::tf::siso::Integrator::UniquePtr integrator_;

public:
  /*===============CONSTRUCTORS===================*/

  //! Constructor
  /*!
      \param Ts sampling time
  */
  PI(double Ts, double Kp = 0.0, double Ki = 0.0) : Ts_(Ts), Kp_(Kp), Ki_(Ki)
  {
    set_parameters();
  }

  void set_parameters()
  {
    if (Ki_ == 0.0)
    {
      integrator_.reset();
    }
    else
    {
      integrator_ = std::make_unique<uclv::systems::tf::siso::Integrator>(Ts_, 1.0);
    }
  }

  //! Destructor
  virtual ~PI() override = default;

  //! Copy Constructor
  PI(const PI& pi) = default;

  virtual PI* clone() const override
  {
    return new PI(this->Ts_, this->Kp_, this->Ki_);  // TODO COPY STATE
  }

  /*==============================================*/

  /*=============GETTER===========================*/
  //! Get the last output
  virtual double getLastOutput() const override
  {
    return last_output_;
  }
  /*==============================================*/

  /*=============SETTER===========================*/
  inline virtual void setTs(double Ts)
  {
    Ts_ = Ts;
    set_parameters();
  }

  //! Change the state so that the output is "output"
  /*!
      \param output output after the state change
  */
  inline virtual void setIntegratorOutput(double output)
  {
    last_output_ = output;
    integrator_->setOutput(output);
  }

  inline virtual void reset()
  {
    last_output_ = 0.0;
    integrator_->reset();
  }

  /*==============================================*/

  /*=============RUNNER===========================*/
  inline virtual double step(double uk) override
  {
    last_output_ = (Kp_ * uk) + (Ki_ * integrator_->step(uk));
    return last_output_;
  }
  /*==============================================*/

  /*=============VARIE===========================*/

  virtual void display() const override
  {
    std::cout << "PI CONTROLLER:\n"
              << "   Ts: " << Ts_ << "\n"
              << "   Kp: " << Kp_ << "\n"
              << "   Ki: " << Ki_ << "\n";
  }

  /*==============================================*/
};

/*=============STATIC FUNS===========================*/
/*==============================================*/

}  // namespace uclv::systems::controllers
