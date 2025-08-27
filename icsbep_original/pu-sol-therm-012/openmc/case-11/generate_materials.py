import openmc

mats = openmc.Materials()

mat = openmc.Material(1)
mat.name = "Plutonium nitrate solution (21.7 g/L)"
mat.set_density('sum')
mat.add_nuclide('Pu239', 4.06270e-05)
mat.add_nuclide('Pu240', 1.02943e-05)
mat.add_nuclide('Pu241', 3.05154e-06)
mat.add_nuclide('Pu242', 6.16420e-07)
mat.add_nuclide('Am241', 3.30688e-07)
mat.add_element('N', 1.41749e-03)
mat.add_element('O', 3.54903e-02)
mat.add_nuclide('H1', 6.36500e-02)
mat.add_element('Fe', 5.39161e-06)
mat.add_element('Cr', 1.73381e-06)
mat.add_element('Ni', 1.22885e-06)
mat.add_s_alpha_beta('c_H_in_H2O')
mats.append(mat)

mat = openmc.Material(2)
mat.name = "Air"
mat.set_density('sum')
mat.add_nuclide('O16', 1.0784e-05)
mat.add_nuclide('N14', 4.3090e-05)
mats.append(mat)

mat = openmc.Material(3)
mat.name = "Stainless steel"
mat.set_density('sum')
mat.add_element('Fe', 6.1344e-02)
mat.add_element('Cr', 1.6472e-02)
mat.add_element('Ni', 8.1050e-03)
mats.append(mat)

mat = openmc.Material(4)
mat.name = "Lucoflex"
mat.set_density('sum')
mat.add_element('C', 2.7365e-02)
mat.add_nuclide('H1', 4.1047e-02)
mat.add_element('Cl', 1.3682e-02)
mats.append(mat)

mat = openmc.Material(5)
mat.name = "Water"
mat.set_density('sum')
mat.add_nuclide('H1', 6.6688e-02)
mat.add_element('O', 3.3344e-02)
mat.add_s_alpha_beta('c_H_in_H2O')
mats.append(mat)

mats.export_to_xml()
