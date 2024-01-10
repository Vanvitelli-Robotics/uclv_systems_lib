#include <uclv_systems_lib/controllers/pi.hpp>

int main()
{
  uclv::systems::controllers::PI pi(0.1, 1.0, 1.0);
  std::cout << "Hello World!" << std::endl;
  return 0;
}
