#include <uclv_systems_lib/controllers/pi.hpp>
#include <uclv_systems_lib/sim/state_space_system_simulator.hpp>
#include <uclv_systems_lib/ss/linear_state_space.hpp>
#include <uclv_systems_lib/continuous_time/continuous_time_linear_state_space.hpp>
#include <uclv_systems_lib/discretization/forward_euler.hpp>
#include <uclv_systems_lib/observers/ekf.hpp>

int main()
{
  const int dim_state = 3;
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
  x0 << 10, 0, 5;

  // define the conitnuous time linear state space system
  auto continuous_time_system_ptr =
      std::make_shared<uclv::systems::ContinuousTimeLinearStateSpace<dim_state, dim_input, dim_output>>(A, B, C, D);
  continuous_time_system_ptr->set_state(x0);
  continuous_time_system_ptr->display();

  // discretized system

  // Create the Forward Euler discretized system
  auto discretized_system =
      std::make_shared<uclv::systems::ForwardEuler<dim_state, dim_input, dim_output>>(continuous_time_system_ptr, 0.1);

  discretized_system->display();

  // extended kalman filter
  Eigen::Matrix<double, dim_state, dim_state> W;
  Eigen::Matrix<double, dim_output, dim_output> V;
  W << 0.1, 0, 0, 0.1, 0, 0, 0.1;
  V << 0.1;

  uclv::systems::ExtendedKalmanFilter<dim_state, dim_input, dim_output> ekf(discretized_system, W, V);

    Eigen::Matrix<double, dim_input, 1> u_k;
    u_k << 1, -1;

    Eigen::Matrix<double, dim_output, 1> y_k;
    y_k = discretized_system->step(u_k);

    Eigen::Matrix<double, dim_state, 1> x_hat_k_k;
    Eigen::Matrix<double, dim_output, 1> y_hat_k;

    ekf.kf_apply(u_k, y_k, W, V, x_hat_k_k, y_hat_k);

    std::cout << "u_k: " << u_k << std::endl;
    std::cout << "y_k: " << y_k << std::endl;

    std::cout << "Estimation result!" << std::endl;
    std::cout << "x_hat_k_k: " << x_hat_k_k << std::endl;
    std::cout << "y_hat_k: " << y_hat_k << std::endl;

  return 0;
}
