// Copyright (c) d-SEAMS developers
//
// Minimal chemfiles mmCIF round trip. chemfiles 0.10.4 crashes in
// mmCIFFormat::read when this file is read back.
#include <chemfiles.hpp>
#include <iostream>

int main() {
  chemfiles::Frame frame;
  frame.set_cell(chemfiles::UnitCell(chemfiles::Vector3D(10, 11, 12)));
  frame.add_atom(chemfiles::Atom("O"), chemfiles::Vector3D(1, 2, 3));
  const char *path = "round.mmcif";
  {
    chemfiles::Trajectory out(path, 'w', "mmCIF");
    out.write(frame);
  }
  chemfiles::Trajectory in(path, 'r', "mmCIF");
  auto got = in.read();
  std::cout << "atoms " << got.size() << "\n";
  return got.size() == 1 ? 0 : 1;
}
