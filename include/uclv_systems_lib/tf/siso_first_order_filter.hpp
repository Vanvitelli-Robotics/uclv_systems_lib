/*
    FirstOrderFilter Class

    First order Low_Pass_Filter

    Copyright 2023 Università della Campania Luigi Vanvitelli

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

/*! \file TF_FIRST_ORDER_FILTER.h
    \brief This class represents a first order filter as Discrete Time Transfer Function System.
*/

#include "siso.hpp"

namespace uclv::systems::tf::siso
{

//!  FirstOrderFilter class: represents a first order filter as Discrete Time Transfer Function System.
/*!
    This class is a first order filter as Discrete Time Transfer Function System.

    It is discretized using the tustin method.

    \sa Linear_System_Interface, TF_INTEGRATOR, TF_SISO, TF_MIMO
*/
class FirstOrderFilter : public SISO
{
private:
  FirstOrderFilter();  // NO DEFAULT CONSTRUCTOR

protected:
  double Ts_;  //!< Sampling time
  double cut_freq_;

  //! Additional filter gain
  double gain_;

public:
  typedef std::shared_ptr<FirstOrderFilter> SharedPtr;
  typedef std::shared_ptr<const FirstOrderFilter> ConstSharedPtr;
  typedef std::weak_ptr<FirstOrderFilter> WeakPtr;
  typedef std::weak_ptr<const FirstOrderFilter> ConstWeakPtr;
  typedef std::unique_ptr<FirstOrderFilter> UniquePtr;

  /*===============CONSTRUCTORS===================*/

  //! Constructor
  /*!
      \param cut_freq cut frequency
      \param Ts sampling time
      \param gain gain of the filter, default = 1
  */
  FirstOrderFilter(double cut_freq, double Ts, double gain = 1.0)
    : SISO(compute_num_coeff(cut_freq, Ts), compute_den_coeff(cut_freq, Ts)), Ts_(Ts), gain_(gain)
  {
    assert(cut_freq > 0.0 && Ts > 0.0 && "FirstOrderFilter: cut_freq and Ts must be > 0.0");
  }

  //! Copy constructor
  FirstOrderFilter(const FirstOrderFilter& tf) = default;

  void set_parameters(double cut_freq, double Ts, double gain = -1.0)
  {
    cut_freq_ = cut_freq < 0.0 ? cut_freq_ : cut_freq;
    Ts_ = Ts;
    gain_ = gain < 0.0 ? gain_ : gain;
    b_vec_ = compute_num_coeff(cut_freq_, Ts);
    a_vec_ = compute_den_coeff(cut_freq_, Ts);
  }

  //! Desctructor
  virtual ~FirstOrderFilter() override = default;

  virtual FirstOrderFilter* clone() const override
  {
    return new FirstOrderFilter(*this);
  }

  /*==============================================*/

  /*===============STATIC FUNCTIONS FOR SIMPLE CONSTRUCTOR WRITING=======*/

  //! INTERNAL, numerator coefficients
  static Eigen::Vector2d compute_num_coeff(double cut_freq, double Ts)
  {
    double tau = 1.0 / (2.0 * M_PI * cut_freq);
    double _alpha = 2.0 * tau / Ts;
    return (1.0 / (1.0 + _alpha)) * Eigen::Vector2d::Ones();
  }

  //! INTERNAL, denominator coefficients
  static Eigen::Vector2d compute_den_coeff(double cut_freq, double Ts)
  {
    double tau = 1.0 / (2.0 * M_PI * cut_freq);
    double _alpha = 2.0 * tau / Ts;
    return Eigen::Vector2d(1.0, ((1.0 - _alpha) / (1.0 + _alpha)));
  }

  /*==============================================*/

  /*=============GETTER===========================*/
  /*==============================================*/

  /*=============SETTER===========================*/

  //! Change the state so that the output is "output"
  /*!
      \param output output after the state change
  */
  inline virtual void setOutput(double output)
  {
    u_vec_[0] = output;
    u_vec_[1] = output;
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
    std::cout << "TF_FIRST_ORDER_FILTER:\n"
              << "   Ts: " << Ts_ << "\n"
              << "   cut_freq: " << (1.0 / (((1.0 / b_vec_[0]) - 1.0) * Ts_ / 2.0)) / (2.0 * M_PI) << "\n"
              << "   cut_freq(dbg): " << cut_freq_ << "\n"
              << "   gain: " << gain_ << "\n";
  }

  /*==============================================*/
};

/*=============STATIC FUNS===========================*/
/*==============================================*/

}  // namespace uclv::systems::tf::siso
