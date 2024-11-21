#include <uclv_systems_lib/controllers/pi.hpp>
#include <uclv_systems_lib/sim/state_space_system_simulator.hpp>
#include <uclv_systems_lib/ss/linear_state_space.hpp>
#include <uclv_systems_lib/continuous_time/continuous_time_linear_state_space.hpp>
#include <uclv_systems_lib/discretization/forward_euler.hpp>

int main()
{
  const int dim_state = 1;
  const int dim_input = 1;
  const int dim_output = 1;

  // define A, B, C, D matrices
  Eigen::Matrix<double, dim_state, dim_state> A;
  Eigen::Matrix<double, dim_state, dim_input> B;
  Eigen::Matrix<double, dim_output, dim_state> C;
  Eigen::Matrix<double, dim_output, dim_input> D;

  // A << -1, 0, 0, 0, -1, 0, 0, 0, -1;
  // B << 1, 0, 0, 1, 0, 0;

  // C << 1, 0, 0;

    // D << 0, 0;
  A << -1;
  B << 0;
  C << 0;
  D << 0;
  Eigen::Matrix<double, dim_state, 1> x0;
  x0 << 1;

  // define the conitnuous time linear state space system
  uclv::systems::ContinuousTimeLinearStateSpace<dim_state, dim_input, dim_output, 1> continuous_time_system(A, B, C, D);
  continuous_time_system.set_state(x0);
  continuous_time_system.display();

  auto discretized_system =
      std::make_shared<uclv::systems::ForwardEuler<dim_state, dim_input, dim_output>>(continuous_time_system, 0.1);
  discretized_system->display();

  uclv::systems::StateSpaceSystemSimulator<dim_state, dim_input, dim_output> simulator(discretized_system);

  Eigen::Matrix<double, dim_input, 1> u_k;
  u_k << 0;

  for (int i = 0; i < 100; i++)
  {
    simulator.simulate(u_k);
  }

  simulator.disp_sim();

  std::cout << "Hello World!" << std::endl;
  return 0;
}
