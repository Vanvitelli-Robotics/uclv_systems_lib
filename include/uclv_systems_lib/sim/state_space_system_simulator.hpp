#pragma once
#include "system_simulator.hpp"
#include "../ss/state_space_interface.hpp"

namespace uclv::systems
{

template <typename Scalar_t, int dim1_state, int dim1_input, int dim1_output, int dim2_state = 1, int dim2_input = 1, int dim2_output = 1>
class StateSpaceSystemSimulator : public SystemSimulator<Scalar_t, dim1_input, dim1_output, dim2_input, dim2_output>
{
public:
  typedef std::shared_ptr<StateSpaceSystemSimulator> SharedPtr;
  typedef std::shared_ptr<const StateSpaceSystemSimulator> ConstSharedPtr;
  typedef std::weak_ptr<StateSpaceSystemSimulator> WeakPtr;
  typedef std::weak_ptr<const StateSpaceSystemSimulator> ConstWeakPtr;
  typedef std::unique_ptr<StateSpaceSystemSimulator> UniquePtr;

  typedef ::uclv::systems::StateSpaceInterface<Scalar_t, dim1_state, dim1_input, dim1_output, dim2_state, dim2_input, dim2_output>
      StateSpaceInterface;

public:
  StateSpaceSystemSimulator() = default;

  StateSpaceSystemSimulator(typename StateSpaceInterface::SharedPtr state_space_system_ptr)
    : SystemSimulator<Scalar_t, dim1_input, dim1_output, dim2_input, dim2_output>(state_space_system_ptr)
    , state_space_system(state_space_system_ptr){};

  virtual ~StateSpaceSystemSimulator() = default;

  virtual inline void reset_history() override
  {
    SystemSimulator<Scalar_t, dim1_input, dim1_output, dim2_input, dim2_output>::reset_history();
    state_history_.clear();
  }

  virtual inline void reserve_history(std::size_t size) override
  {
    SystemSimulator<Scalar_t, dim1_input, dim1_output, dim2_input, dim2_output>::reserve_history(size);
    state_history_.reserve(size);
  }

  const std::vector<Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& get_state_history() const
  {
    return state_history_;
  }

  void reset(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>>& x0)
  {
    reset_history();
    state_space_system->set_state(x0);
  }

  const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>& get_initial_state() const
  {
    if (state_history_.empty())
    {
      return state_space_system->get_state();
    }
    return x0_;
  }

  const Eigen::Matrix<Scalar_t, dim1_state, dim2_state>& get_final_state() const
  {
    if (state_history_.empty())
    {
      return state_space_system->get_state();
    }
    return state_history_.back();
  }

  virtual inline void simulate(const Eigen::Ref<const Eigen::Matrix<Scalar_t, dim1_input, dim2_input>>& u_k) override
  {
    assert(system != nullptr && "[SystemSimulator] System not set");

    // if first step save initial condition
    if (state_history_.empty())
    {
      x0_ = state_space_system->get_state();
    }

    SystemSimulator<Scalar_t, dim1_input, dim1_output, dim2_input, dim2_output>::simulate(u_k);
    state_history_.push_back(state_space_system->get_state());
  }

  virtual void disp_sim() const override
  {
    SystemSimulator<Scalar_t, dim1_input, dim1_output, dim2_input, dim2_output>::disp_sim();
    std::cout << "State History: \n";
    for (std::size_t i = 0; i < state_history_.size(); i++)
    {
      std::cout << "x[" << i << "]:\n" << state_history_[i] << "\n";
    }
  }

public:
  typename StateSpaceInterface::SharedPtr state_space_system;

protected:
  std::vector<Eigen::Matrix<Scalar_t, dim1_state, dim2_state>> state_history_;
  Eigen::Matrix<Scalar_t, dim1_state, dim2_state> x0_;
};
}  // namespace uclv::systems
