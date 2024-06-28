#include <uclv_systems_lib/controllers/pi.hpp>
#include <uclv_systems_lib/sim/state_space_system_simulator.hpp>
#include <uclv_systems_lib/ss/linear_state_space.hpp>

int main()
{
  uclv::systems::controllers::PI pi(0.1, 1.0, 1.0);
  uclv::systems::StateSpaceSystemSimulator<2, 3, 4, 5, 6, 7> simulator1;
  auto system = std::make_shared<uclv::systems::LinearStateSpace<2, 3, 4, 5>>();
  uclv::systems::StateSpaceSystemSimulator<2, 3, 4, 5, 5, 5> simulator2(system);
  std::cout << "Hello World!" << std::endl;
  return 0;
}
