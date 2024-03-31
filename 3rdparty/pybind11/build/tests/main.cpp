#include <filesystem>
int main(int argc, char ** argv) {
  std::filesystem::path p(argv[0]);
  return p.string().length();
}