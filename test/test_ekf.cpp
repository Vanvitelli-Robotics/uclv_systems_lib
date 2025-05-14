#include <uclv_systems_lib/controllers/pi.hpp>
#include <uclv_systems_lib/sim/state_space_system_simulator.hpp>
#include <uclv_systems_lib/ss/linear_state_space.hpp>
#include <uclv_systems_lib/continuous_time/continuous_time_linear_state_space.hpp>
#include <uclv_systems_lib/discretization/forward_euler.hpp>
#include <uclv_systems_lib/observers/ekf.hpp>
#define SCALAR_TYPE double

int main()
{
  const int dim_state = 2;
  const int dim_input = 2;
  const int dim_output = 2;

  // define A, B, C, D matrices
  Eigen::Matrix<double, dim_state, dim_state> A;
  Eigen::Matrix<double, dim_state, dim_input> B;
  Eigen::Matrix<double, dim_output, dim_state> C;
  Eigen::Matrix<double, dim_output, dim_input> D;

  A << -1, 0, 0, -1;
  B << 1, 0, 0, 1;

  C << 1, 0, 0, 2;

  D << 0, 0;
  // A << -1;
  // B << 0;
  // C << 0;
  // D << 0;
  Eigen::Matrix<double, dim_state, 1> x0;
  x0 << 10, 0;

  // define the conitnuous time linear state space system
  auto continuous_time_system_ptr =
      std::make_shared<uclv::systems::ContinuousTimeLinearStateSpace<SCALAR_TYPE,dim_state, dim_input, dim_output>>(A, B, C, D,x0);

  // Create the Forward Euler discretized system
  auto discretized_system =
      std::make_shared<uclv::systems::ForwardEuler<SCALAR_TYPE,dim_state, dim_input, dim_output>>(continuous_time_system_ptr, 0.1);
  discretized_system->set_state(x0);
  discretized_system->display();

  // extended kalman filter
  // define W and V covariance matrices
  Eigen::Matrix<double, dim_state, dim_state> W;
  Eigen::Matrix<double, dim_output, dim_output> V;
  W << 0.1, 0, 0, 0.1;
  V << 0.1, 0, 0, 0.1;

  // start from a different initial state
  Eigen::Matrix<double, dim_state, 1> x0_hat;
  x0_hat << -0, 0;

  uclv::systems::ExtendedKalmanFilter<SCALAR_TYPE,dim_state, dim_input, dim_output> ekf(discretized_system, W, V);
  ekf.set_state(x0_hat);

  // simulation of the system and the observer
  Eigen::Matrix<double, dim_input, 1> u_k;
  u_k << 1, -1;
  Eigen::Matrix<double, dim_output, 1> y_k;
  Eigen::Matrix<double, dim_state, 1> x_hat_k_k;
  Eigen::Matrix<double, dim_output, 1> y_hat_k;
  uclv::systems::StateSpaceSystemSimulator<SCALAR_TYPE,dim_state, dim_input, dim_output> simulator(discretized_system);

  for (int i = 0; i < 20; i++)
  {
    simulator.simulate(u_k);
    y_k = discretized_system->get_output();
    ekf.kf_apply(u_k, y_k, W, V);
    x_hat_k_k = ekf.get_state();
    y_hat_k = ekf.get_output();
    std::cout << "Real state!\n" << std::endl;

    std::cout << "x_k: " << discretized_system->get_state().transpose() << std::endl;
    std::cout << "y_k: " << y_k.transpose() << std::endl;
    std::cout << "Estimation result!\n" << std::endl;
    std::cout << "x_hat_k_k: " << x_hat_k_k.transpose() << std::endl;
    std::cout << "y_hat_k: " << y_hat_k.transpose() << std::endl;
  }

  ekf.display();

  return 0;
}
