#include "HODLR_Matrix.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

int add(int i, int j) 
{
    return i + j;
}

Eigen::MatrixXd inv(const Eigen::MatrixXd &xs)
{
  return xs.inverse();
}

class Pet
{
    public:
        Pet(const std::string &name, int hunger) : name(name), hunger(hunger) {}
        ~Pet() {}

        void go_for_a_walk() { hunger++; }
        const std::string &get_name() const { return name; }
        int get_hunger() const { return hunger; }

    private:
        std::string name;
        int hunger;
};

class Animal 
{
public:
    virtual ~Animal() { }
    virtual std::string go(int n_times) = 0;
};

class PyAnimal : public Animal 
{
public:
    /* Inherit the constructors */
    using Animal::Animal;

    /* Trampoline (need one for each virtual function) */
    std::string go(int n_times) override {
        PYBIND11_OVERLOAD_PURE(
            std::string, /* Return type */
            Animal,      /* Parent class */
            go,          /* Name of function in C++ (must match Python name) */
            n_times      /* Argument(s) */
        );
    }
};

class Dog : public Animal 
{
    public:
    std::string go(int n_times) override 
    {
        std::string result;
        for (int i=0; i<n_times; ++i)
            result += "woof! ";
        return result;
    }
};

std::string call_go(Animal *animal) {
    return animal->go(3);
}

PYBIND11_MODULE(pyhodlrlib, m) 
{
    m.doc() = "This is the Python Wrapper to the HODLRlib Kernel Object";

    // define add function
    m.def("add", &add, "A function which adds two numbers");

    m.def("inv", &inv);

    py::class_<Animal, PyAnimal /* <--- trampoline*/> animal(m, "Animal");
    animal
        .def(py::init<>())
        .def("go", &Animal::go);

    // bindings to Pet class
    py::class_<Pet>(m, "Pet")
        .def(py::init<const std::string &, int>())
        .def("go_for_a_walk", &Pet::go_for_a_walk)
        .def("get_hunger", &Pet::get_hunger)
        .def("get_name", &Pet::get_name);

    m.def("call_go", &call_go);
}
