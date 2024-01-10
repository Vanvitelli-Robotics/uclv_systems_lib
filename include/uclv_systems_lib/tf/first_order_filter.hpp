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

#include "siso_first_order_filter.hpp"
#include <vector>

namespace uclv::systems::tf::mimo
{

//!  FirstOrderFilter class: represents a first order filter as Discrete Time Transfer Function System.
/*!
    This class is a first order filter as Discrete Time Transfer Function System.

    It is discretized using the tustin method.

    \sa Linear_System_Interface, TF_INTEGRATOR, TF_SISO, TF_MIMO
*/
template <int FILTER_SIZE = Eigen::Dynamic>
class FirstOrderFilter
{
private:
  FirstOrderFilter();  // NO DEFAULT CONSTRUCTOR

protected:
  double Ts_;  //!< Sampling time
  double cut_freq_;

  //! Additional filter gain
  double gain_;

  std::vector<uclv::systems::tf::siso::FirstOrderFilter::SharedPtr> filters_;
  Eigen::Matrix<double, FILTER_SIZE, 1> output_;

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
  FirstOrderFilter(unsigned int size, double cut_freq, double Ts, double gain = 1.0) : Ts_(Ts), gain_(gain)
  {
    if (FILTER_SIZE != Eigen::Dynamic)
    {
      assert(size == FILTER_SIZE &&
             "FirstOrderFilter::FirstOrderFilter: cannot set a fixed size different from the template size");
    }
    set_parameters(cut_freq, Ts, size, gain);
  }

  FirstOrderFilter(double cut_freq, double Ts, double gain = 1.0) : Ts_(Ts), gain_(gain)
  {
    if (FILTER_SIZE == Eigen::Dynamic)
    {
      throw std::runtime_error("FirstOrderFilter::FirstOrderFilter: you have to specify the size of the filter");
    }
    set_parameters(cut_freq, Ts, gain, FILTER_SIZE);
  }

  FirstOrderFilter(unsigned int size) : FirstOrderFilter(size, 1.0, 1.0, 0.0)
  {
  }

  //! Copy constructor
  FirstOrderFilter(const FirstOrderFilter& tf) = default;

  void set_parameters(double cut_freq = -1, double Ts = -1, double gain = -1.0, int size = -1)
  {
    size = size < 0 ? filters_.size() : size;
    if (FILTER_SIZE != Eigen::Dynamic)
    {
      assert(size == FILTER_SIZE &&
             "FirstOrderFilter::FirstOrderFilter: cannot set a fixed size different from the template size");
    }

    cut_freq_ = cut_freq < 0 ? cut_freq_ : cut_freq;
    Ts_ = Ts < 0 ? Ts_ : Ts;
    gain_ = gain < 0 ? gain_ : gain;
    output_.resize(size);

    filters_.reserve(size);

    if ((size_t)size < filters_.size())
    {
      filters_.resize(size);
    }

    for (auto& filter : filters_)
    {
      filter->set_parameters(cut_freq_, Ts_, gain_);
    }

    for (unsigned int i = filters_.size(); i < (size_t)size; i++)
    {
      filters_.push_back(std::make_shared<siso::FirstOrderFilter>(cut_freq_, Ts_, gain_));
    }
  }

  //! Desctructor
  virtual ~FirstOrderFilter() = default;

  virtual FirstOrderFilter* clone() const
  {
    return new FirstOrderFilter(*this);
  }

  /*==============================================*/

  /*===============STATIC FUNCTIONS FOR SIMPLE CONSTRUCTOR WRITING=======*/

  /*==============================================*/

  /*=============GETTER===========================*/
  /*==============================================*/

  /*=============SETTER===========================*/

  //! Change the state so that the output is "output"
  /*!
      \param output output after the state change
  */
  inline virtual void setOutput(const Eigen::Matrix<double, FILTER_SIZE, 1>& output)
  {
    assert(output.size() == (Eigen::Index)filters_.size() && "FirstOrderFilter::setOutput: output size mismatch");

    for (size_t i = 0; i < filters_.size(); i++)
    {
      filters_[i]->setOutput(output[i]);
    }
    output_ = output;
  }
  /*==============================================*/

  /*=============RUNNER===========================*/
  inline virtual const Eigen::Matrix<double, FILTER_SIZE, 1>& step(const Eigen::Matrix<double, FILTER_SIZE, 1>& uk)
  {
    assert(uk.size() == (Eigen::Index)filters_.size() && "FirstOrderFilter::step: input size mismatch");

    for (size_t i = 0; i < filters_.size(); i++)
    {
      output_[i] = filters_[i]->step(uk[i]);
    }

    return output_;
  }
  /*==============================================*/

  /*=============VARIE===========================*/
  inline virtual void reset()
  {
    for (auto& filter : filters_)
    {
      filter->reset();
    }
    output_.setZero();
  }

  virtual unsigned int getSizeInput() const
  {
    return filters_.size();
  }

  virtual unsigned int getSizeOutput() const
  {
    return filters_.size();
  }

  //! Get the last output
  virtual const Eigen::Matrix<double, FILTER_SIZE, 1>& getLastOutput() const
  {
    return output_;
  }

  virtual void display() const
  {
    display_tf();
  }

  virtual void display_tf() const
  {
    std::cout << "TF_FIRST_ORDER_FILTER:\n";
    std::cout << "------------------------------------\n";
    for (size_t i = 0; i < filters_.size(); i++)
    {
      std::cout << "Filter " << i << ":\n";
      filters_[i]->display_tf();
      std::cout << "------------------------------------\n";
    }
    std::cout << "####################################\n";
  }

  /*==============================================*/
};

/*=============STATIC FUNS===========================*/
/*==============================================*/

}  // namespace uclv::systems::tf::mimo
