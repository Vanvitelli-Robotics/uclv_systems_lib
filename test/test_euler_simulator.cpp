#include <uclv_systems_lib/controllers/pi.hpp>
#include <uclv_systems_lib/sim/state_space_system_simulator.hpp>
#include <uclv_systems_lib/ss/linear_state_space.hpp>
#include <uclv_systems_lib/continuous_time/continuous_time_linear_state_space.hpp>
#include <uclv_systems_lib/discretization/forward_euler.hpp>

int main()
{const int dim_state = 3;
  const int dim_input = 2;
  const int dim_output = 1;

  // define A, B, C, D matrices
  Eigen::Matrix<double, dim_state, dim_state> A;
  Eigen::Matrix<double, dim_state, dim_input> B;
  Eigen::Matrix<double, dim_output, dim_state> C;
  Eigen::Matrix<double, dim_output, dim_input> D;

  A << -1, 0, 0, 0, -1, 0, 0, 0, -1;
  B << 1, 0, 0, 1, 0, 0;

  C << 1, 0, 0;

  D << 0, 0;
  // A << -1;
  // B << 0;
  // C << 0;
  // D << 0;
  Eigen::Matrix<double, dim_state, 1> x0;
  x0 << 1, 0, 0;

  // define the conitnuous time linear state space system
  auto continuous_time_system_ptr =
      std::make_shared<uclv::systems::ContinuousTimeLinearStateSpace<dim_state, dim_input, dim_output>>(A, B, C, D);
  continuous_time_system_ptr->display();

  // discretized system

  // Create the Forward Euler discretized system
  auto discretized_system =
      std::make_shared<uclv::systems::ForwardEuler<dim_state, dim_input, dim_output>>(continuous_time_system_ptr, 0.1);
  discretized_system->set_state(x0);
  discretized_system->display();

  Eigen::Matrix<double, dim_input, 1> u_k;
  u_k << 0, -0;


  uclv::systems::StateSpaceSystemSimulator<dim_state, dim_input, dim_output> simulator(discretized_system);

  for (int i = 0; i < 10; i++)
  {
    simulator.simulate(u_k);
  }

  simulator.disp_sim();

  std::cout << "Hello World!" << std::endl;

  simulator.system->display();


  // test jacobian functions
  Eigen::Matrix<double, dim_state, dim_state> jacobx;
  Eigen::Matrix<double, dim_output, dim_state> jacobx_output;
  Eigen::Matrix<double, dim_output, dim_input> jacobu_output;
  Eigen::Matrix<double, dim_state, dim_input> jacobu;

  continuous_time_system_ptr->jacobx_state_fcn(continuous_time_system_ptr->get_state(), u_k, jacobx);
  continuous_time_system_ptr->jacobu_state_fcn(continuous_time_system_ptr->get_state(), u_k, jacobu);
  continuous_time_system_ptr->jacobx_output_fcn(continuous_time_system_ptr->get_state(), u_k, jacobx_output);
  continuous_time_system_ptr->jacobu_output_fcn(continuous_time_system_ptr->get_state(), u_k, jacobu_output);

  std::cout << "jacobx: \n" << jacobx << std::endl;
  std::cout << "jacobu: \n" << jacobu << std::endl;
  std::cout << "jacobx_output: \n" << jacobx_output << std::endl;
  std::cout << "jacobu_output: \n" << jacobu_output << std::endl;

  discretized_system->jacobx_state_fcn(discretized_system->get_state(), u_k, jacobx);
  discretized_system->jacobu_state_fcn(discretized_system->get_state(), u_k, jacobu);
  discretized_system->jacobx_output_fcn(discretized_system->get_state(), u_k, jacobx_output);
  discretized_system->jacobu_output_fcn(discretized_system->get_state(), u_k, jacobu_output);

  std::cout << "jacobx: \n" << jacobx << std::endl;
  std::cout << "jacobu: \n" << jacobu << std::endl;
  std::cout << "jacobx_output: \n" << jacobx_output << std::endl;
  std::cout << "jacobu_output: \n" << jacobu_output << std::endl;

  return 0;
}
