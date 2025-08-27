import openmc

mats = openmc.Materials()

mat = openmc.Material(1)
mat.name = "Plutonium solution"
mat.set_density('sum')
mat.add_nuclide('Pu238', 3.3429e-07)
mat.add_nuclide('Pu239', 6.7826e-05)
mat.add_nuclide('Pu240', 7.0032e-05)
mat.add_nuclide('Pu241', 1.7699e-05)
mat.add_nuclide('Pu242', 7.6195e-06)
mat.add_nuclide('Am241', 1.7570e-06)
mat.add_element('Gd', 2.2195e-08)
mat.add_nuclide('H1', 6.1242e-02)
mat.add_element('N', 2.1108e-03)
mat.add_element('O', 3.6227e-02)
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
