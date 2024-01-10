/*
    TF_SISO Class Discrete Time Transfer Function SISO (Linear Filter)

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

#include "../siso_interface.hpp"

#include <Eigen/Dense>
#include <memory>
#include <iostream>

/*! \file siso.hpp
    \brief This class represents a generic Discrete Time Transfer Function System.
*/

/*
         b0 + b1*z-1 + b2*z-2 + ... + bn*z-n
    yk = ------------------------------------ uk
         a0 + a1*z-1 + a2*z-2 + ... + am*z-m

    yk = ([b0 b1 b2 ... bn]/a0)*[uk uk-1 uk-2 ... uk-n]T - ([a1 a2 ... am]/a0)*[yk-1 yk-2 ... yk-m]T

    b_vec = [b0 b1 b2 ... bn]T / a0  -> size n+1

    a_vec = [a1 a2 ... am]T / a0 -> size m

    u_vec = [uk uk-1 uk-2 ... uk-n]T

    y_vec = [yk-1 yk-2 ... yk-m]^T

    yk = b_vecT * u_vec - a_vecT * y_vec

*/

namespace uclv::systems::tf::siso
{
//!  TF_SISO class: represents a generic Discrete Time Transfer Function System.
/*!
    This class is a generic Discrete Time Transfer Function System in the form:

    \verbatim

    Z transform:

            b0 + b1*z^-1 + b2*z^-2 + ... + bn*z^-n
    Y(z) = ---------------------------------------- U(z)
            a0 + a1*z^-1 + a2*z^-2 + ... + am*z^-m

    Discrete Time (time k,  ^T = transpose):

    y(k) = ([b0 b1 b2 ... bn]/a0)*[u(k) u(k-1) u(k-2) ... u(k-n)]^T - ([a1 a2 ... am]/a0)*[y(k-1) y(k-2) ... y(k-m)]^T

    Short formulation:

    b_vec = [b0 b1 b2 ... bn]^T / a0  -> size n+1

    a_vec = [a1 a2 ... am]^T / a0 -> size m

    u_vec = [u(k) u(k-1) u(k-2) ... u(k-n)]^T

    y_vec = [y(k-1) y(k-2) ... y(k-m)]^T

    yk = b_vec^T * u_vec - a_vec^T * y_vec

    \endverbatim

    It stores the internal system state

    \sa Linear_System_Interface, TF_INTEGRATOR, TF_FIRST_ORDER_FILTER, TF_MIMO
*/
class SISO : public SISOInterface
{
public:
  typedef std::shared_ptr<SISO> SharedPtr;
  typedef std::shared_ptr<const SISO> ConstSharedPtr;
  typedef std::weak_ptr<SISO> WeakPtr;
  typedef std::weak_ptr<const SISO> ConstWeakPtr;
  typedef std::unique_ptr<SISO> UniquePtr;

  // STATIC

  //! INTERNAL simplify numerator
  /*!
    This function transforms the numerator vector [b0 b1 ... bn 0 0 ... 0]T in b_vec=[b0 b1 ... bn]T / a0
    \param num_coeff numerator coefficients
    \param den_coeff denominator coefficients
    \return simplified numerator coefficients
  */
  inline static Eigen::VectorXd simplifyNumerator(const Eigen::VectorXd& num_coeff, const Eigen::VectorXd& den_coeff)
  {
    if (den_coeff[0] == 0)
    {
      throw std::domain_error("[tf::SISO::simplifyNumerator] Division by 0");
    }

    if (num_coeff.size() == 0)  // case void vector
    {
      Eigen::VectorXd semplified_numerator(1);
      semplified_numerator.setZero();
      return semplified_numerator;
    }

    // find first nonzero coeff
    unsigned int effective_size = num_coeff.size();
    while (effective_size != 1)
    {
      if (num_coeff[effective_size - 1] != 0.0)
        break;
      effective_size--;
    }
    return (num_coeff.segment(0, effective_size) / den_coeff[0]);
  }

  //
  //! INTERNAL simplify and reduce the denominator
  /*!
    This function transforms the denominator [a0 a1 ... am 0 0 ... 0]T in a_vec=[a1 a2 ... am]T / a0
    \param den_coeff denominator coefficients
    \return simplified denominator coefficients
  */
  inline static Eigen::VectorXd simplifyAndReduceDenominator(const Eigen::VectorXd& den_coeff)
  {
    if (den_coeff[0] == 0)
    {
      throw std::domain_error("[tf::SISO::simplifyAndReduceDenominator] Division by 0");
    }

    // find first nonzero coeff
    unsigned int effective_size = den_coeff.size();
    while (effective_size != 1)
    {
      if (den_coeff[effective_size - 1] != 0.0)
        break;
      effective_size--;
    }

    if (effective_size == 1)
    {
      Eigen::VectorXd semplified_denominator(0);
      return semplified_denominator;  // return a void vector
    }

    return (den_coeff.segment(1, effective_size - 1) / den_coeff[0]);
  }

private:
protected:
  //! Numerator coefficients, size n+1
  Eigen::VectorXd b_vec_;  // n+1
  //! Reduced Denominator coefficients, without a0, size m
  Eigen::VectorXd a_vec_;  // m

  //! State vector, last inputs, size n+1
  Eigen::VectorXd u_vec_;  // n+1
  //! State vector, last outputs, size m
  Eigen::VectorXd y_vec_;  // m
  //! very last output, size 1
  double y_k_ = 0;  // last output

public:
  /*===============CONSTRUCTORS===================*/

  //! Contructor
  /*!
    Construct the discrete time transfer function

    \param num_coeff numerator coefficients [b0, ...., bn]
    \param den_coeff denominator coefficients [a0, ...., an]
  */
  SISO(const Eigen::VectorXd& num_coeff, const Eigen::VectorXd& den_coeff)
    : b_vec_(simplifyNumerator(num_coeff, den_coeff)), a_vec_(simplifyAndReduceDenominator(den_coeff))
  {
    u_vec_.resize(b_vec_.size());
    y_vec_.resize(a_vec_.size());

    u_vec_.setZero();
    y_vec_.setZero();
  }

  //! Zero Contructor
  /*!
    Construct a ZERO discrete time transfer function
  */
  SISO() : SISO(Eigen::VectorXd{ { 0.0 } }, Eigen::VectorXd{ { 1.0 } })
  {
  }

  //! Copy Constructor
  SISO(const SISO& tf) = default;

  virtual ~SISO() = default;

  //! Clone the object
  virtual SISO* clone() const
  {
    return new SISO(*this);
  }

  /*==============================================*/

  /*=============GETTER===========================*/

  //! Get the numerator order n
  inline virtual unsigned int getNumeratorOrder() const
  {
    return (b_vec_.size() - 1);
  }

  //! Get the denominator order m
  inline virtual unsigned int getDenominatorOrder() const
  {
    return (a_vec_.size());
  }

  /*==============================================*/

  /*=============SETTER===========================*/

  /*==============================================*/

  /*=============RUNNER===========================*/

  inline virtual double step(double u_k)
  {
    // Shift u_vec
    const unsigned int no = u_vec_.size() - 1;
    if (no != 0)
    {
      Eigen::VectorXd tmp_u_vec = u_vec_;
      u_vec_.segment(1, no) = tmp_u_vec.segment(0, no);
    }
    u_vec_[0] = u_k;

    // Shift y_vec
    const unsigned int deno = y_vec_.size();
    if (deno == 0)
    {
      y_k_ = b_vec_.dot(u_vec_);
      return y_k_;
    }
    if (deno != 1)
    {
      Eigen::VectorXd tmp_y_vec = y_vec_;
      y_vec_.segment(1, deno - 1) = tmp_y_vec.segment(0, deno - 1);
    }
    y_vec_[0] = y_k_;

    y_k_ = (b_vec_.dot(u_vec_)) - (a_vec_.dot(y_vec_));

    return y_k_;
  }

  /*==============================================*/

  /*=============VARIE===========================*/
  inline virtual void reset()
  {
    u_vec_.setZero();
    if (y_vec_.size() != 0)
      y_vec_.setZero();
    y_k_ = 0.0;
  }

  virtual unsigned int getSizeInput() const
  {
    return 1;
  }

  virtual unsigned int getSizeOutput() const
  {
    return 1;
  }

  //! Get the last output
  virtual double getLastOutput() const
  {
    return y_k_;
  }

  virtual void display() const
  {
    display_tf();
  }

  //! Display the transfer function on the std out
  virtual void display_tf() const
  {
    std::cout << "TF_SISO:\n"
              << "   Num_Coeff= " << b_vec_ << "\n"
              << "   Den_Coeff= 1.0 " << a_vec_ << "\n"
              << "   State"
              << "\n"
              << "   u_vec= " << u_vec_ << "\n"
              << "   y_vec= " << y_vec_ << "\n"
              << "   y_k= " << y_k_ << "\n";
  }

  /*==============================================*/
};

}  // namespace uclv::systems::tf::siso
