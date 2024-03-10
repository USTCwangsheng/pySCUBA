py::class_<geo::SASAParameter>(m, "SASAParameter")
        .def(py::init<>(), py::arg("solventradius") = 1.4, py::arg("soluteradius") = 1.4,py::arg("nsurfacepoints")=256)
        .def_readwrite("solventradius", &geo::SASAParameter::solventradius)
        .def_readwrite("soluteradius",&geo::SASAParameter::soluteradius)
        .def_readwrite("nsurfacepoints",&geo::SASAParameter::nsurfacepoints);
py::class_<geo::AtomSASA>(m, "AtomSASA")
        .def(py::init<>())
        .def(py::init<const geo::SASAParameter &>)
        .def(py::init<const geo::LocalFrame &, const geo::SASAParameter &>())
        .def_readwrite("solutecrd", &geo::AtomSASA::solutecrd_)
        .def_readwrite("surfacepoints",&geo::AtomSASA::surfacepoints)
        .def("updateAtomSASA",&geo::updateAtomSASA);