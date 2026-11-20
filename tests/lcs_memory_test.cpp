#include <cplib/str/lcs.hpp>
#include <cstdlib>
#include <new>
#include <string>

namespace {
bool measureAllocations = false;
std::size_t allocatedBytes = 0;
}
void *operator new(std::size_t size) {
    if (measureAllocations)
        allocatedBytes += size;
    if (void *memory = std::malloc(size == 0 ? 1 : size))
        return memory;
    throw std::bad_alloc();
}
void operator delete(void *memory) noexcept {
    std::free(memory);
}
void operator delete(void *memory, std::size_t) noexcept {
    std::free(memory);
}
int main() {
    constexpr std::size_t n = 4000;
    std::string a(n, 'a'), b(n, 'a');
    allocatedBytes = 0;
    measureAllocations = true;
    auto length = cplib::LCS(a, b);
    measureAllocations = false;
    assert(length == cplib::Int(n));
    assert(allocatedBytes <= 4 * (n + 1) * sizeof(cplib::Int));
    assert(cplib::LCS(std::string{}, b) == 0);
    assert(cplib::LCS(std::vector<int>{1, 2, 1, 3}, std::vector<int>{2, 1, 4, 3}) == 3);
}
