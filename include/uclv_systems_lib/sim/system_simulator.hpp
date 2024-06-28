#pragma once
#include "../system_interface.hpp"
#include <vector>

namespace uclv::systems
{

template <int dim1_input, int dim1_output, int dim2_input = 1, int dim2_output = 1>
class SystemSimulator
{
public:
  typedef std::shared_ptr<SystemSimulator> SharedPtr;
  typedef std::shared_ptr<const SystemSimulator> ConstSharedPtr;
  typedef std::weak_ptr<SystemSimulator> WeakPtr;
  typedef std::weak_ptr<const SystemSimulator> ConstWeakPtr;
  typedef std::unique_ptr<SystemSimulator> UniquePtr;

  typedef ::uclv::systems::SystemInterface<dim1_input, dim1_output, dim2_input, dim2_output> SystemInterface;

protected:
  std::vector<Eigen::Matrix<double, dim1_input, dim2_input>> input_history_;
  std::vector<Eigen::Matrix<double, dim1_output, dim2_output>> output_history_;

public:
  typename SystemInterface::SharedPtr system;

  SystemSimulator(){};
  SystemSimulator(typename SystemInterface::SharedPtr system_ptr) : system(system_ptr){};

  virtual ~SystemSimulator() = default;

  virtual inline void reset_history()
  {
    input_history_.clear();
    output_history_.clear();
  }

  virtual inline void reserve_history(std::size_t size)
  {
    input_history_.reserve(size);
    output_history_.reserve(size);
  }

  const std::vector<Eigen::Matrix<double, dim1_input, dim2_input>>& get_input_history() const
  {
    return input_history_;
  }

  const std::vector<Eigen::Matrix<double, dim1_output, dim2_output>>& get_output_history() const
  {
    return output_history_;
  }

  virtual inline void simulate(const Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>& u_k)
  {
    assert(system != nullptr && "[SystemSimulator] System not set");
    input_history_.push_back(u_k);
    output_history_.push_back(system->step(u_k));
  }

  inline void simulate(const std::vector<Eigen::Matrix<double, dim1_input, dim2_input>>& input_vector)
  {
    for (const auto& u_k : input_vector)
    {
      simulate(u_k);
    }
  }

  inline void simulate(const std::vector<Eigen::Ref<const Eigen::Matrix<double, dim1_input, dim2_input>>>& input_vector)
  {
    for (const auto& u_k : input_vector)
    {
      simulate(u_k);
    }
  }

  inline void simulate(const std::vector<Eigen::Ref<Eigen::Matrix<double, dim1_input, dim2_input>>>& input_vector)
  {
    for (const auto& u_k : input_vector)
    {
      simulate(u_k);
    }
  }

  virtual void disp_sim() const
  {
    std::cout << "System Simulation Report: \n";
    std::cout << "Input History: \n";
    for (std::size_t i = 0; i < input_history_.size(); i++)
    {
      std::cout << "u[" << i << "]:\n" << input_history_[i] << "\n";
    }
    std::cout << "Output History: \n";
    for (std::size_t i = 0; i < output_history_.size(); i++)
    {
      std::cout << "y[" << i << "]:\n" << output_history_[i] << "\n";
    }
  }
};
}  // namespace uclv::systems
