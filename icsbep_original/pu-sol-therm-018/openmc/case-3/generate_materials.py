import openmc

mats = openmc.Materials()

mat = openmc.Material(1)
mat.name = "Plutonium solution"
mat.set_density('sum')
mat.add_nuclide('Pu238', 5.0954e-07)
mat.add_nuclide('Pu239', 1.0338e-04)
mat.add_nuclide('Pu240', 1.0675e-04)
mat.add_nuclide('Pu241', 2.6978e-05)
mat.add_nuclide('Pu242', 1.1614e-05)
mat.add_nuclide('Am241', 2.6781e-06)
mat.add_element('Gd', 3.3832e-08)
mat.add_nuclide('H1', 5.8391e-02)
mat.add_element('N', 3.2152e-03)
mat.add_element('O', 3.7736e-02)
mat.add_s_alpha_beta('c_H_in_H2O')
mats.append(mat)

mat = openmc.Material(2)
mat.name = "steel clad"
mat.set_density('sum')
mat.add_element('Fe', 6.3278e-02)
mat.add_element('Cr', 1.6532e-02)
mat.add_element('Ni', 6.5095e-03)
mats.append(mat)

mat = openmc.Material(3)
mat.name = "water"
mat.set_density('sum')
mat.add_nuclide('H1', 6.6691e-02)
mat.add_element('O', 3.3346e-02)
mat.add_s_alpha_beta('c_H_in_H2O')
mats.append(mat)

mats.export_to_xml()
