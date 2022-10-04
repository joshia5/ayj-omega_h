#include <Omega_h_library.hpp>
#include <Omega_h_irrule.hpp>

using namespace Omega_h;

void test_irrule(Library *lib) {
  auto comm = lib->world();

  auto irrule = TetGaussLobatto(4);

  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  test_irrule(&lib);

  return 0;
}
